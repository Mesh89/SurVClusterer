#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include <chrono>

#include "libs/IntervalTree.h"
#include "libs/kseq.h"
#include "libs/ssw_cpp.h"
#include "libs/cxxopts.h"

#include "common.h"

std::unordered_map<std::string, std::string> chrs;
int max_prec_dist, max_imprec_dist;
double min_prec_frac_overlap, min_imprec_frac_overlap;
int max_prec_len_diff, max_imprec_len_diff;

StripedSmithWaterman::Aligner aligner(1,4,6,1);
StripedSmithWaterman::Filter filter;

bool ignore_seq = false;

int distance(sv_t& sv1, sv_t& sv2) {
	if (sv1.chr != sv2.chr) return INT32_MAX;
	if ((sv1.type == "DUP" && sv2.type == "INS") || (sv1.type == "INS" && sv2.type == "DUP")) {
		return std::min(abs(sv1.start-sv2.start), abs(sv1.end-sv2.end));
	}
	return std::max(abs(sv1.start-sv2.start), abs(sv1.end-sv2.end));
}
double overlap(const sv_t& sv1, const sv_t& sv2) {
	if (sv1.type == "INS" && sv2.type == "INS") return 1.0;
	if (sv1.end == sv1.start || sv2.end == sv2.start) return 1.0;
    int overlap_bp = std::max(0, std::min(sv1.end, sv2.end)-std::max(sv1.start, sv2.start));
    return overlap_bp/double(std::min(sv1.end-sv1.start, sv2.end-sv2.start));
}
int len_diff(sv_t& sv1, sv_t& sv2) {
	if (sv1.type == "DUP" && sv2.type == "DUP") return abs(sv1.len()-sv2.len());
	else if (sv1.type == "INS" && sv2.type == "INS") return abs((int) (sv1.ins_seq.length()-sv2.ins_seq.length()));
	return 0;
}

bool is_compatible_dup_dup(sv_t& sv1, sv_t& sv2) {
	if (!sv1.precise || !sv2.precise) {
		return  distance(sv1, sv2) <= max_imprec_dist &&
				overlap(sv1, sv2) >= min_imprec_frac_overlap &&
				abs(sv1.len()-sv2.len()) <= max_imprec_len_diff;
	} else {
		return  distance(sv1, sv2) <= max_prec_dist &&
				overlap(sv1, sv2) >= min_prec_frac_overlap &&
				abs(sv1.len()-sv2.len()) <= max_prec_len_diff;
	}
}

bool check_ins_dup_seq(sv_t& ins_sv, sv_t& dup_sv) {
	if (ignore_seq) return true;

	int max_len_diff = (ins_sv.precise && dup_sv.precise) ? max_prec_len_diff : max_imprec_len_diff;
	if (dup_sv.len() > ins_sv.ins_seq.length()+max_len_diff) return false;

	std::string dup_seq = chrs[dup_sv.chr].substr(dup_sv.start, dup_sv.end-dup_sv.start);
	std::string ext_dup_seq;
	while (ext_dup_seq.length() <= ins_sv.ins_seq.length()) ext_dup_seq += dup_seq;

	if (ext_dup_seq.length() > UINT16_MAX || ins_sv.ins_seq.length() > UINT16_MAX) return true; // TODO: ssw score is a uint16_t, so we cannot compare strings longer than that

	StripedSmithWaterman::Alignment alignment;
	aligner.Align(ins_sv.ins_seq.data(), ext_dup_seq.data(), ext_dup_seq.length(), filter, &alignment);
	return alignment.query_end-alignment.query_begin >= ins_sv.ins_seq.length()*0.8;
}

bool is_compatible_ins_dup(sv_t& sv1, sv_t& sv2) {
	if (!sv1.precise || !sv2.precise) {
		return distance(sv1, sv2) <= max_imprec_dist && check_ins_dup_seq(sv1, sv2);
	} else {
		return distance(sv1, sv2) <= max_prec_dist && check_ins_dup_seq(sv1, sv2);
	}
}

bool check_ins_ins_seq(sv_t& sv1, sv_t& sv2) {
	if (ignore_seq) return true;
	if (sv1.ins_seq.length() > UINT16_MAX || sv2.ins_seq.length() > UINT16_MAX) return true; // TODO: ssw score is a uint16_t, so we cannot compare strings longer than that

	int max_len_diff = (sv1.precise && sv2.precise) ? max_prec_len_diff : max_imprec_len_diff;
	if (abs((int) (sv1.ins_seq.length()-sv2.ins_seq.length())) > max_len_diff) return false;

	sv_t& l_sv = sv1.start < sv2.start ? sv1 : sv2;
	sv_t& r_sv = sv1.start < sv2.start ? sv2 : sv1;
	std::string extra_seq = chrs[sv1.chr].substr(l_sv.start, r_sv.start-l_sv.start);

	std::string lsv_seq = l_sv.ins_seq + extra_seq;
	std::string rsv_seq = extra_seq + r_sv.ins_seq;

	StripedSmithWaterman::Alignment alignment;
	std::string& query = lsv_seq.length() < rsv_seq.length() ? lsv_seq : rsv_seq;
	std::string& ref = lsv_seq.length() < rsv_seq.length() ? rsv_seq : lsv_seq;
	aligner.Align(query.data(), ref.data(), ref.length(), filter, &alignment);
	return alignment.query_end-alignment.query_begin >= query.length()*0.8;
}

bool is_compatible_ins_ins(sv_t& sv1, sv_t& sv2) {
	if (!sv1.precise || !sv2.precise) {
		return distance(sv1, sv2) <= max_imprec_dist && check_ins_ins_seq(sv1, sv2);
	} else {
		return distance(sv1, sv2) <= max_prec_dist && check_ins_ins_seq(sv1, sv2);
	}
}

bool check_ins_seq(sv_t& sv1, sv_t& sv2) {
	if (sv1.type == "DUP" && sv2.type == "DUP") {
		return true;
	} else if (sv1.type == "DUP" && sv2.type == "INS") {
		return check_ins_dup_seq(sv2, sv1);
	} else if (sv1.type == "INS" && sv2.type == "DUP") {
		return check_ins_dup_seq(sv1, sv2);
	} else {
		return check_ins_ins_seq(sv1, sv2);
	}
}

bool is_compatible(sv_t& sv1, sv_t& sv2) {
	if (sv1.type == "DUP" && sv2.type == "DUP") {
		return is_compatible_dup_dup(sv1, sv2);
	} else if (sv1.type == "DUP" && sv2.type == "INS") {
		return is_compatible_ins_dup(sv2, sv1);
	} else if (sv1.type == "INS" && sv2.type == "DUP") {
		return is_compatible_ins_dup(sv1, sv2);
	} else {
		return is_compatible_ins_ins(sv1, sv2);
	}
}

int main(int argc, char* argv[]) {

	cxxopts::Options options("compare-ins", "Given a VCF or SV file with benchmark insertions and one with the called ones, reports for "
											"each benchmark insertion whether it was called.");

	options.add_options()
		("benchmark_file", "VCF or SV file with the benchmark calls.", cxxopts::value<std::string>())
		("called_file", "VCF or SV file with the calls to benchmark.", cxxopts::value<std::string>())
		("T,tandem-repeats", "Tandem repeat file in BED format.", cxxopts::value<std::string>())
		("R,reference", "Reference file in fasta format.", cxxopts::value<std::string>())
		("d,max_dist_precise", "Maximum distance allowed between the breakpoints of two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("D,max_dist_imprecise", "Maximum distance allowed between the breakpoints when at least one of the variants is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("v,min_overlap_precise", "Minimum overlap required between two precise variants, expressed as the fraction of the shortest variant"
				" that is covered by the longest variant.", cxxopts::value<double>()->default_value("0.8"))
		("V,min_overlap_imprecise", "Minimum overlap required between two variants when at least one is imprecise, "
				"expressed as the fraction of the shortest variant that is covered by the longest variant.",
				cxxopts::value<double>()->default_value("0.5"))
		("s,max_len_diff_precise", "Maximum length difference allowed between two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("S,max_len_diff_imprecise", "Maximum length difference allowed when at least one variant is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("r,report", "Print report only.", cxxopts::value<bool>()->default_value("false"))
		("f,print-fps", "Print false positive IDs to stderr.", cxxopts::value<bool>()->default_value("false"))
		("I,ignore-seq", "Only compare insertions by position.", cxxopts::value<bool>()->default_value("false"))
		("i,force-ids", "Generate new IDs for the variants, required when variants do not have IDs or have duplicated IDs.",
				cxxopts::value<bool>()->default_value("false"))
		("a,all-imprecise", "Treat all deletions as imprecise.", cxxopts::value<bool>()->default_value("false"))
		("h,help", "Print usage");

	options.parse_positional({"benchmark_file", "called_file"});
	options.positional_help("benchmark_file called_file");
	options.show_positional_help();
	options.allow_unrecognised_options();
	auto parsed_args = options.parse(argc, argv);

	if (parsed_args.count("help") || !parsed_args.count("benchmark_file") || !parsed_args.count("called_file")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

	std::string benchmark_fname = parsed_args["benchmark_file"].as<std::string>();
	std::vector<sv_t> benchmark_svs = read_sv_list(benchmark_fname.c_str());
	std::string called_fname = parsed_args["called_file"].as<std::string>();
	std::vector<sv_t> called_svs = read_sv_list(called_fname.c_str());
	max_prec_dist = parsed_args["max_dist_precise"].as<int>();
	max_imprec_dist = parsed_args["max_dist_imprecise"].as<int>();
	min_prec_frac_overlap = parsed_args["min_overlap_precise"].as<double>();
	min_imprec_frac_overlap = parsed_args["min_overlap_imprecise"].as<double>();
	max_prec_len_diff = parsed_args["max_len_diff_precise"].as<int>();
	max_imprec_len_diff = parsed_args["max_len_diff_imprecise"].as<int>();
	bool report = parsed_args["report"].as<bool>(), print_fp = parsed_args["print-fps"].as<bool>();
	ignore_seq = parsed_args["ignore-seq"].as<bool>();
	if (parsed_args["force-ids"].as<bool>()) {
		for (int j = 0; j < benchmark_svs.size(); j++) benchmark_svs[j].id = "INS_" + std::to_string(j);
		for (int j = 0; j < called_svs.size(); j++) called_svs[j].id = "INS_" + std::to_string(j);
	}
	if (parsed_args["all-imprecise"].as<bool>()) {
		max_prec_dist = max_imprec_dist;
		min_prec_frac_overlap = min_imprec_frac_overlap;
		max_prec_len_diff = max_imprec_len_diff;
	}

	if (!parsed_args.count("reference") && !ignore_seq) {
		std::cerr << "If no reference is provided, the --ignore-seq flag must be provided." << std::endl;
		exit(0);
	}

	auto is_ins_func = [](sv_t& sv) {return (sv.type != "INS" && sv.type != "DUP") || sv.len() <= 0;};
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_ins_func), benchmark_svs.end());
	called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_ins_func), called_svs.end());

	std::unordered_map<std::string, std::vector<sv_t> > called_svs_by_chr;
	for (sv_t& sv : called_svs) {
		if (!ignore_seq && !sv.has_seq() && sv.type == "INS") continue;
		called_svs_by_chr[sv.chr].push_back(sv);
	}

	if (parsed_args.count("reference")) {
		FILE* fastaf = fopen(parsed_args["reference"].as<std::string>().c_str(), "r");
		kseq_t *seq = kseq_init(fileno(fastaf));
		while (kseq_read(seq) >= 0) {
			chrs[std::string(seq->name.s)] = seq->seq.s;
		}
		kseq_destroy(seq);
	}

	std::unordered_map<std::string, IntervalTree<repeat_t>*> reps_i;
    if (parsed_args.count("tandem-repeats")) {
		std::unordered_map<std::string, std::vector<repeat_t>> reps;
		std::unordered_map<std::string, std::vector<Interval<repeat_t>>> reps_iv;
		std::string line;
		std::ifstream rep_f(parsed_args["tandem-repeats"].as<std::string>());
		while (getline(rep_f, line)) {
			if (line[0] == '#') continue;
			repeat_t r(line);
			reps[r.chr].push_back(r);
			reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
		}

		for (auto it = reps.begin(); it != reps.end(); it++) {
			reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
		}
    }

    /* == Compare benchmark svs with called svs == */
    std::set<std::string> b_tps, c_tps;

    for (sv_t& bsv : benchmark_svs) {
        std::vector<repeat_t> reps_containing_bsv;
		if (reps_i[bsv.chr] != NULL) {
			std::vector<Interval<repeat_t>> intervals_temp = reps_i[bsv.chr]->findOverlapping(bsv.start, bsv.end);
			for (auto& iv : intervals_temp) {
				repeat_t rep = iv.value;
				if (rep.contains(bsv)) {
					reps_containing_bsv.push_back(rep);
				}
			}
		}

		bool matched = false;
        for (sv_t& csv : called_svs_by_chr[bsv.chr]) {
            if (is_compatible(bsv, csv)) {
                if (!report) std::cout << bsv.id << " " << csv.id << std::endl;
                b_tps.insert(bsv.id);
                c_tps.insert(csv.id);
                matched = true;
            }
            if (matched) continue;

            int max_len_diff = (bsv.precise && csv.precise) ? max_prec_len_diff : max_imprec_len_diff;
			if (len_diff(bsv, csv) <= max_len_diff) {
				bool same_tr = false;
				for (repeat_t& rep : reps_containing_bsv) {
					if (rep.intersects(csv)) {
						same_tr = true;
						break;
					}
				}

				if (same_tr && check_ins_seq(bsv, csv)) {
					if (!report) std::cout << bsv.id << " " << csv.id << " REP" << std::endl;
					b_tps.insert(bsv.id);
					c_tps.insert(csv.id);
					matched = true;
				}
			}
        }


        if (!matched) {
        	if (!report) std::cout << bsv.id << " NONE" << std::endl;
        }
	}

	if (report) {
		std::cout.precision(2);
		std::cout << "RECALL: " << b_tps.size() << "/" << benchmark_svs.size() << " = " << double(b_tps.size())/benchmark_svs.size() << std::endl;
		std::cout << "PRECISION: " << c_tps.size() << "/" << called_svs.size() << " = " << double(c_tps.size())/called_svs.size() << std::endl;
	}

	if (print_fp)
	for (sv_t& sv : called_svs) {
		if (!c_tps.count(sv.id)) {
			std::cerr << sv.id << std::endl;
		}
	}
}
