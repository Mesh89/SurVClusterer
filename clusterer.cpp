#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unistd.h>

#include "htslib/vcf.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read)


struct chr_seq_t {
    char* seq;
    hts_pos_t len;

    chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
    ~chr_seq_t() {delete[] seq;}
};
struct chr_seqs_map_t {
    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname) {
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strcpy(chr_seq, seq->seq.s);
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};


std::unordered_map<std::string, int> sample2id;
std::vector<std::string> sample_names;
chr_seqs_map_t chr_seqs;


std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_ins_seq(bcf_hdr_t* hdr, bcf1_t* sv) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		return sv->d.allele[1];
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) return data;

	return "";
}

struct sv_t {
    std::string id, chr, type, ins_seq, sample;
    int start, end;
    char str1, str2;
    bool precise;

    sv_t(std::string id, std::string chr, int start, char str1, int end, char str2, std::string type, std::string seq,
         std::string sample, std::string precise)
    : id(id), chr(chr), start(start), str1(str1), end(end), str2(str2), type(type), ins_seq(seq), sample(sample),
    precise(precise == "PRECISE") {}

    sv_t(bcf_hdr_t* hdr, bcf1_t* vcf_sv, std::string& sample) : sample(sample), str1('N'), str2('N') {
    	bcf_unpack(vcf_sv, BCF_UN_INFO);
    	id = vcf_sv->d.id;
    	chr = bcf_seqname(hdr, vcf_sv);
    	start = vcf_sv->pos;
    	end = get_sv_end(hdr, vcf_sv);
    	type = get_sv_type(hdr, vcf_sv);
    	ins_seq = get_ins_seq(hdr, vcf_sv);
    	precise = true;
    	if (bcf_get_info_flag(hdr, vcf_sv, "IMPRECISE", NULL, 0) == 1) precise = false;
    }

    int len() {
        if (type == "INS") return ins_seq.length();
        else return end-start;
    }

    std::string to_string() {
        std::stringstream ss;
        ss << sample << "." << id << "," << chr << ":" << start << "-" << end << "," << len() << ",";
        ss << (precise ? "PRECISE" : "IMPRECISE");
        return ss.str();
    }

    std::string to_coordinates() {
    	std::stringstream ss;
    	ss << chr << ":" << start << "-" << end;
    	return ss.str();
    }
};
bool operator < (const sv_t& sv1, const sv_t& sv2) {
    if (sv1.chr != sv2.chr) return sv1.chr < sv2.chr;
    return sv1.start < sv2.start;
}
int distance(const sv_t& sv1, const sv_t& sv2) {
    if (sv1.chr != sv2.chr) return INT32_MAX;
    return abs(sv1.start-sv2.start) + abs(sv1.end-sv2.end);
}

struct idx_size_t {
    int idx, size;

    idx_size_t() : idx(-1), size(0) {}
    idx_size_t(int idx, int size) : idx(idx), size(size) {}
};
bool operator > (const idx_size_t& s1, const idx_size_t& s2) {
    return s1.size > s2.size;
};

std::vector<sv_t> svs;
std::vector<std::vector<int>> compatibility_list;
int max_dist, min_overlap;

int overlap(const sv_t& sv1, const sv_t& sv2) {
    return std::max(0, std::min(sv1.end, sv2.end)-std::max(sv1.start, sv2.start));
}

bool is_compatible(const sv_t& sv1, const sv_t& sv2) {
    return distance(sv1, sv2) <= max_dist && overlap(sv1, sv2) >= min_overlap && abs((int)(sv1.ins_seq.length()-sv2.ins_seq.length())) <= max_dist;
}

int median(std::vector<int>& v) {
    if (v.size()%2 == 0) {
        return (v[v.size()/2-1] + v[v.size()/2])/2;
    } else {
        return v[v.size()/2];
    }
}

std::vector<std::vector<int>> compute_minimal_clique_cover(int start, int end) {
    std::vector<idx_size_t> idx_by_neighborhood_size;
    for (int i = start; i < end; i++) {
        idx_by_neighborhood_size.emplace_back(i, compatibility_list[i].size());
    }
    std::sort(idx_by_neighborhood_size.begin(), idx_by_neighborhood_size.end(), std::greater<idx_size_t>());

    std::vector<std::vector<int>> cliques;
    int cliques_counter = 0;
    std::vector<int> cliques_idx(svs.size(), -1);
    for (idx_size_t& curr_idx : idx_by_neighborhood_size) {
        // find cliqued neighbors
        std::vector<idx_size_t> cliqued_neighbors;
        for (int adj_idx : compatibility_list[curr_idx.idx]) {
            int clique_idx = cliques_idx[adj_idx];
            if (clique_idx != -1 && is_compatible(svs[curr_idx.idx], svs[adj_idx])) {
                cliqued_neighbors.emplace_back(adj_idx, cliques[clique_idx].size());
            }
        }

        // check if I can join any neighboring cliques
        std::sort(cliqued_neighbors.begin(), cliqued_neighbors.end(), std::greater<idx_size_t>());
        bool used = false;
        for (idx_size_t& neighbor_idx : cliqued_neighbors) {
            int neighboor_clique_idx = cliques_idx[neighbor_idx.idx];
            std::vector<int> &neighboor_clique = cliques[neighboor_clique_idx];
            bool belongs_to_clique = true;
            for (int clique_elem_idx : neighboor_clique) {
                if (!is_compatible(svs[curr_idx.idx], svs[clique_elem_idx])) {
                    belongs_to_clique = false;
                    break;
                }
            }

            if (belongs_to_clique) {
                neighboor_clique.push_back(curr_idx.idx);
                used = true;
                break;
            }
        }

        if (!used) {
            cliques.push_back(std::vector<int>());
            cliques[cliques_counter].push_back(curr_idx.idx);
            cliques_idx[curr_idx.idx] = cliques_counter;
            cliques_counter++;
        }
    }

    return cliques;
}

sv_t choose_sv(std::vector<sv_t>& svs) {
    std::unordered_map<std::string, int> counts;
    for (sv_t& sv : svs) {
        if (sv.precise) {
            counts[sv.to_coordinates()]++;
        }
    }
    if (counts.empty()) { // if no precise accepts imprecises
        for (sv_t& sv : svs) {
            counts[sv.to_coordinates()]++;
        }
    }

    int max_freq = 0;
    std::string str;
    for (auto& e : counts) {
        if (e.second > max_freq) {
            max_freq = e.second;
            str = e.first;
        }
    }

    for (sv_t& sv : svs) {
        if (sv.to_coordinates() == str) return sv;
    }
    std::cerr << "SHOULD NOT BE HERE." << std::endl;
    exit(1);
}

int cluster_id = 0;
void print_cliques(std::vector<std::vector<int>>& cliques, htsFile* out_vcf_file, bcf_hdr_t* out_hdr, std::ofstream& out_sv_file, chr_seqs_map_t& chr_seqs) {
	int n_samples = bcf_hdr_nsamples(out_hdr);
	const char** coos = new const char*[n_samples];
	const char** sizes = new const char*[n_samples]; // could be int, but I did not find an easy way to set a variable number of ints in the format
	const char** prec = new const char*[n_samples];

	std::vector<std::pair<sv_t, int> > sv_order;
	for (auto& clique : cliques) {
		std::vector<sv_t> clique_svs;
		for (int idx : clique) {
			clique_svs.push_back(svs[idx]);
		}
		sv_t chosen_sv = choose_sv(clique_svs);
		sv_order.push_back({chosen_sv, sv_order.size()});
	}
	std::sort(sv_order.begin(), sv_order.end(),
			[out_hdr](const std::pair<sv_t, int>& sv1, const std::pair<sv_t, int>& sv2) {
		int cid1 = bcf_hdr_name2id(out_hdr, sv1.first.chr.c_str());
		int cid2 = bcf_hdr_name2id(out_hdr, sv2.first.chr.c_str());
		return std::tie(cid1, sv1.first.start) < std::tie(cid2, sv2.first.start);
	});

	bcf1_t* vcf_sv = bcf_init();
	for (auto& clique_sv_idx : sv_order) {
		int clique_idx = clique_sv_idx.second;
		std::vector<int>& clique = cliques[clique_idx];

        std::set<std::string> unique_samples;
        std::vector<sv_t> clique_svs;
        for (int idx : clique) {
            unique_samples.insert(svs[idx].sample);
            clique_svs.push_back(svs[idx]);
        }

        sv_t& chosen_sv = clique_sv_idx.first;

        bcf_clear(vcf_sv);

        // set basic info
        vcf_sv->rid = bcf_hdr_name2id(out_hdr, chosen_sv.chr.c_str());
        vcf_sv->pos = chosen_sv.start;
        std::string id = "CLUSTER_" + std::to_string(cluster_id);
        bcf_update_id(out_hdr, vcf_sv, id.c_str());
        char* chr_seq = chr_seqs.get_seq(chosen_sv.chr);
        std::string alleles = std::string(1, chr_seq[chosen_sv.start]) + ",<" + chosen_sv.type + ">";
		bcf_update_alleles_str(out_hdr, vcf_sv, alleles.c_str());

		// set INFO
		int int_conv = chosen_sv.end+1;
		bcf_update_info_int32(out_hdr, vcf_sv, "END", &int_conv, 1);
		bcf_update_info_string(out_hdr, vcf_sv, "SVTYPE", chosen_sv.type.c_str());
		int_conv = unique_samples.size();
		bcf_update_info_int32(out_hdr, vcf_sv, "N_SAMPLES", &int_conv, 1);
		int_conv = clique.size();
		bcf_update_info_int32(out_hdr, vcf_sv, "N_SVS", &int_conv, 1);

		// set GT
		int* gts = new int[n_samples];
		std::fill(gts, gts+n_samples, bcf_gt_unphased(0));
		for (sv_t& sv : clique_svs) {
			gts[sample2id[sv.sample]] = bcf_gt_unphased(1);
		}
		bcf_update_genotypes(out_hdr, vcf_sv, gts, n_samples);

		// set CO (original coordinates)
		std::vector<std::string> coos_str(n_samples);
		std::vector<std::string> sizes_str(n_samples);
		std::vector<std::string> prec_str(n_samples);
		for (sv_t& sv : clique_svs) {
			if (!coos_str[sample2id[sv.sample]].empty()) {
				coos_str[sample2id[sv.sample]] += ",";
				sizes_str[sample2id[sv.sample]] += ",";
				prec_str[sample2id[sv.sample]] += ",";
			}
			coos_str[sample2id[sv.sample]] += std::to_string(sv.start) + "-" + std::to_string(sv.end);
			sizes_str[sample2id[sv.sample]] += std::to_string(sv.end-sv.start);
			prec_str[sample2id[sv.sample]] += (sv.precise ? "P" : "I");
		}
		for (int i = 0; i < n_samples; i++) {
			if (coos_str[i].empty()) {
				coos_str[i] = ".";
				sizes_str[i] = ".";
				prec_str[i] = ".";
			}
			coos[i] = coos_str[i].c_str();
			sizes[i] = sizes_str[i].c_str();
			prec[i] = prec_str[i].c_str();
		}
		bcf_update_format_string(out_hdr, vcf_sv, "CO", coos, n_samples);
		bcf_update_format_string(out_hdr, vcf_sv, "SZ", sizes, n_samples);
		bcf_update_format_string(out_hdr, vcf_sv, "IP", prec, n_samples);

		if (bcf_write(out_vcf_file, out_hdr, vcf_sv) != 0) {
        	throw std::runtime_error("Failed to write record to " + std::string(out_vcf_file->fn));
        }

        out_sv_file << "CLUSTER_" << cluster_id << " " << chosen_sv.chr << " " << chosen_sv.start << " " << chosen_sv.str1;
        out_sv_file << " " << chosen_sv.chr << " " << chosen_sv.end << " " << chosen_sv.str2 << " " << chosen_sv.type << " ";
        out_sv_file << unique_samples.size() << " " << clique.size() << " " << (chosen_sv.ins_seq.empty() ? "NA" : chosen_sv.ins_seq) << " ";

        cluster_id++;

        std::sort(clique_svs.begin(), clique_svs.end());
        for (sv_t sv : clique_svs) {
        	out_sv_file << " " << sv.to_string();
        }
        out_sv_file << std::endl;
    }
    bcf_destroy(vcf_sv);

    delete[] coos;
    delete[] sizes;
    delete[] prec;
}

bcf_hrec_t* generate_contig_hrec() {
	bcf_hrec_t* contig_hrec = new bcf_hrec_t;
	contig_hrec->type = BCF_HL_CTG;
	contig_hrec->key = strdup("contig");
	contig_hrec->value = NULL;
	contig_hrec->keys = contig_hrec->vals = NULL;
	contig_hrec->nkeys = 0;
	int r1 = bcf_hrec_add_key(contig_hrec, "ID", 2);
	int r2 = bcf_hrec_add_key(contig_hrec, "length", 6);
	if (r1 || r2) {
		throw std::runtime_error("Failed to create contig to VCF header.");
	}
	return contig_hrec;
}
bcf_hdr_t* generate_vcf_header() {

	bcf_hdr_t* out_hdr = bcf_hdr_init("w");

	// add contigs
	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(chr_seqs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(out_hdr, hrec);
	}

	// add INFO tags
	int len;
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the indel (DEL or DUP).\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, end_tag, &len));

	const char* n_samples_tag = "##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description=\"Number of samples having at least one SV in this cluster.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, n_samples_tag, &len));

	const char* n_svs_tag = "##INFO=<ID=N_SVS,Number=1,Type=Integer,Description=\"Number of SVs clustered. This may be different from N_SAMPLES as more than "
			"one SV per sample may be clustered.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, n_svs_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, gt_tag, &len));

	// add FORMAT tags
	const char* og_sv_tag = "##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Coordinates of the original SV(s).\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, og_sv_tag, &len));

	const char* sz_sv_tag = "##FORMAT=<ID=SZ,Number=1,Type=String,Description=\"Size(s) of the original SV(s).\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, sz_sv_tag, &len));

	const char* pr_sv_tag = "##FORMAT=<ID=IP,Number=1,Type=String,Description=\"Whether the original SV(s) were precise (P) or imprecise (I). "
			"Anything without an explicit IMPRECISE tag is considered precise.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, pr_sv_tag, &len));

	int i = 0;
	for (std::string& sample_name : sample_names) {
		sample2id[sample_name] = i++; // bcf_hdr_nsamples(out_hdr) does not work unless we sync the header
		bcf_hdr_add_sample(out_hdr, sample_name.c_str());
	}

	return out_hdr;
}

int main(int argc, char* argv[]) {

	if (argc < 6) {
		std::cout << "Cluster SVs." << std::endl;
		std::cout << "exec filelist reference max_dist min_overlap out_prefix" << std::endl;
		return 0;
	}

    std::ifstream file_list(argv[1]);
    std::string reference_fname = argv[2];
    max_dist = std::stoi(argv[3]);
    min_overlap = std::stoi(argv[4]);
    std::string out_prefix = argv[5];
    std::string out_vcf_fname = out_prefix + ".vcf.gz";
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    std::string out_sv_fname = out_prefix + ".sv";
    std::ofstream out_sv_file(out_sv_fname);

    chr_seqs.read_fasta_into_map(reference_fname);

    std::string sample_name, sample_sv_fpath;
    bcf1_t* vcf_record = bcf_init();
    while (file_list >> sample_name >> sample_sv_fpath) {
    	std::cout << "Reading SVs from " << sample_sv_fpath << std::endl;
    	htsFile* sample_sv_file = bcf_open(sample_sv_fpath.c_str(), "r");
    	bcf_hdr_t* vcf_header = bcf_hdr_read(sample_sv_file);
		if (vcf_header == NULL) {
			throw std::runtime_error("Failed to read the VCF header of " + sample_sv_fpath + ".");
		}
		// TODO: using header of first file for output. Probably we should implement a consistency check
    	while (bcf_read(sample_sv_file, vcf_header, vcf_record) == 0) {
    		sv_t sv(vcf_header, vcf_record, sample_name);
    		svs.push_back(sv);
    		compatibility_list.push_back(std::vector<int>());
    	}
		bcf_hdr_destroy(vcf_header);
    	hts_close(sample_sv_file);
    	sample_names.push_back(sample_name);
    }
    bcf_destroy(vcf_record);
    std::cout << "Finished reading SVs." << std::endl;

    bcf_hdr_t* out_hdr = generate_vcf_header();
    if (bcf_hdr_write(out_vcf_file, out_hdr) != 0) {
    	throw std::runtime_error("Could not write header to " + std::string(out_vcf_file->fn));
    }

    std::cout << "Clustering " << svs.size() << " SVs." << std::endl;
    std::sort(svs.begin(), svs.end());
    int i_prev = 0;
    for (int i = 0; i < svs.size(); i++) {
        sv_t& sv1 = svs[i];
        for (int j = i+1; j < svs.size(); j++) {
            sv_t& sv2 = svs[j];
            if (sv1.chr != sv2.chr || sv2.start-sv1.start > max_dist) break;
            if (is_compatible(sv1, sv2)) {
                compatibility_list[i].push_back(j);
                compatibility_list[j].push_back(i);
            }
        }
        compatibility_list[i].shrink_to_fit();

        if (i > 0 && (sv1.chr != svs[i-1].chr || sv1.start-svs[i-1].start > max_dist)) {
            std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, i);
            print_cliques(cliques, out_vcf_file, out_hdr, out_sv_file, chr_seqs);

            for (int j = i_prev; j < i; j++) {
                compatibility_list[j].clear();
                compatibility_list[j].shrink_to_fit();
            }
            i_prev = i;
        }

        if (i > 0 && i % 1000000 == 0) std::cout << "Clustered " << i << " SVs." << std::endl;
    }

    std::cout << "Generating minimum clique cover." << std::endl;
    std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, svs.size());

    std::cout << "Writing output." << std::endl;
    print_cliques(cliques, out_vcf_file, out_hdr, out_sv_file, chr_seqs);

    bcf_hdr_destroy(out_hdr);
    hts_close(out_vcf_file);

    std::cout << "Finished." << std::endl;
}
