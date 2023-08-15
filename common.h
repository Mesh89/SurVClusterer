#ifndef SVCOMPARE_COMMON_H
#define SVCOMPARE_COMMON_H

#include <unistd.h>
#include "htslib/vcf.h"
#include "libs/kseq.h"
KSEQ_INIT(int, read)

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

int get_sv_len(bcf_hdr_t* hdr, bcf1_t* sv) {
	int* data = NULL;
	int size = 0;
	bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
	if (size > 0) {
		int len = data[0];
		delete[] data;
		return len;
	}

	if (sv->d.allele[1][0] != '<') {
		return strlen(sv->d.allele[1])-strlen(sv->d.allele[0]);
	}

	return 0;
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

    if (get_sv_type(hdr, sv) == "INS") return sv->pos;

    int svlen = get_sv_len(hdr, sv);
	if (svlen != 0) {
		return sv->pos + abs(svlen);
	}

    std::cerr << "Warning: SV (ID=" << sv->d.id << ") has no END or SVLEN annotation. Assuming END == POS." << std::endl;
    return sv->pos;
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


bool ends_with(const char* str, const char* suffix) {
	int l_string = strlen(str);
	int l_suffix = strlen(suffix);
	return strcmp(str+(l_string-l_suffix), suffix) == 0;
}

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

struct sv_t {
    std::string id, chr, type, ins_seq, sample;
    int start, end;
    int _len = 0;
    bool precise, incomplete_ass;
    int* gts = NULL, n_gts;

    sv_t() : start(0), end(0), precise(true), incomplete_ass(true), n_gts(0) {}
    sv_t(bcf_hdr_t* hdr, bcf1_t* vcf_sv, std::string& sample) : sample(sample) {
    	bcf_unpack(vcf_sv, BCF_UN_INFO);
    	id = vcf_sv->d.id;
    	chr = bcf_seqname(hdr, vcf_sv);
    	start = vcf_sv->pos;
    	end = get_sv_end(hdr, vcf_sv);
    	type = get_sv_type(hdr, vcf_sv);
    	ins_seq = get_ins_seq(hdr, vcf_sv);
    	_len = get_sv_len(hdr, vcf_sv);
    	int n = 0;
    	n_gts = bcf_get_genotypes(hdr, vcf_sv, &gts, &n);
    	precise = true;
    	incomplete_ass = false;
    	if (bcf_get_info_flag(hdr, vcf_sv, "IMPRECISE", NULL, 0) == 1) precise = false;
    	if (bcf_get_info_flag(hdr, vcf_sv, "INCOMPLETE_ASSEMBLY", NULL, 0) == 1) incomplete_ass = true;
    }

    int len() {
    	if (!_len) {
    		if (type == "DEL") return start-end;
    		else return end-start;
    	}
        return _len;
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

    bool has_seq() {
    	return !ins_seq.empty() && ins_seq != "NA";
    }
};
bool operator < (const sv_t& sv1, const sv_t& sv2) {
    if (sv1.chr != sv2.chr) return sv1.chr < sv2.chr;
    return std::tie(sv1.start, sv1.end) < std::tie(sv2.start, sv2.end);
}

std::vector<sv_t> read_sv_list(const char* filename) {
	std::vector<sv_t> svs;
	if (ends_with(filename, ".vcf.gz") || ends_with(filename, ".vcf") || ends_with(filename, ".bcf")) { // vcf format
		htsFile* file = bcf_open(filename, "r");
		bcf_hdr_t* hdr = bcf_hdr_read(file);
		bcf1_t* line = bcf_init();
		while (bcf_read(file, hdr, line) == 0) {
			std::string sample;
			svs.push_back(sv_t(hdr, line, sample));
		}
		bcf_destroy(line);
		bcf_hdr_destroy(hdr);
		bcf_close(file);
	} else {
		std::ifstream file(filename);
		std::string line;
		while (getline(file, line)) {
			std::stringstream ss(line);
			sv_t sv;
			char strand;
			ss >> sv.id >> sv.chr >> sv.start >> strand >> sv.chr >> sv.end >> strand >> sv.type >> sv.ins_seq;
			svs.push_back(sv);
		}
	}
	return svs;
}

struct repeat_t {
    std::string chr;
    int start, end;

    repeat_t() : start(0), end(0) {}
    repeat_t(std::string& line) {
        std::stringstream ss(line);
        ss >> chr >> start >> end;
    }

    bool contains(sv_t& sv) {
        return sv.chr == chr && start <= sv.start && end >= sv.end;
    }
    bool intersects(sv_t& sv) {
        return sv.chr == chr && sv.start <= end && sv.end >= start;
    }
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

#endif //SVCOMPARE_COMMON_H
