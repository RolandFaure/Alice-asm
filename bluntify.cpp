/**
 * @file bluntify.cpp
 * @brief Implements functions for processing and "bluntifying" GFA (Graphical Fragment Assembly) files by removing overlaps, trimming contigs, and cleaning up zero-length contigs.
 *
 * This file provides several utilities for manipulating GFA files, including:
 * - Removing overlaps from links between contigs.
 * - Splitting contigs at overlap breakpoints.
 * - Removing contigs of zero length and updating links accordingly.
 * - Supporting both basic and advanced ("fancier") overlap removal strategies.
 *
 * Main Functions:
 * - reverse_complement: Computes the reverse complement of a DNA sequence.
 * - basic_overlap_removal: Removes overlaps from GFA links and trims contigs accordingly.
 * - fancier_overlap_removal: Splits contigs at overlap breakpoints and updates links, with options for minimum contig length and overlap reporting.
 * - remove_contigs_of_length_0: Removes contigs with zero length and rewires links to maintain graph connectivity.
 * - bluntify: Orchestrates the bluntification process by chaining the above steps and managing intermediate files.
 * - bluntify_main: Command-line interface for the bluntify tool, handling arguments and invoking the main workflow.
 *
 * Internal Data Structures:
 * - Link, ContigSides, LinkKey, Neighbor, WrittenLink: Structures for representing GFA links, contig sides, and neighbor relationships.
 * - Various hash functors for use in unordered containers.
 *
 * Utility Functions:
 * - split_tab, join_tail_with_tabs, split_name_parts: String manipulation helpers for parsing GFA lines.
 * - parse_overlap1, parse_overlap2: Parse CIGAR strings to determine overlap lengths.
 * - orient_from_end1, orient_from_end2: Determine orientation symbols from link ends.
 * - parse_gfa_for_overlap_steps: Parses a GFA file and populates data structures for contigs and links.
 *
 * Usage:
 *   bluntify <input.gfa> <output.gfa> [-n|--no_overlaps] [--tmpdir TMPDIR] [--version]
 *
 * Dependencies:
 *   Requires C++17 or later for <filesystem> and other standard library features.
 *
 * Author: Roland Faure
 * Version: 0.1
 */
#include "bluntify.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using std::getline;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::size_t;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace {

const string kVersion = "0.1";

struct Link {
	string name1;
	int end1 = 0;
	int overlap1 = 0;
	string name2;
	int end2 = 0;
	int overlap2 = 0;
};

struct ContigSides {
	vector<int> left;
	vector<int> right;
};

struct LinkKey {
	string name1;
	int end1 = 0;
	string name2;
	int end2 = 0;

	bool operator==(const LinkKey& other) const {
		return name1 == other.name1 && end1 == other.end1 && name2 == other.name2 && end2 == other.end2;
	}
};

struct LinkKeyHash {
	size_t operator()(const LinkKey& k) const {
		size_t h1 = std::hash<string>{}(k.name1);
		size_t h2 = std::hash<int>{}(k.end1);
		size_t h3 = std::hash<string>{}(k.name2);
		size_t h4 = std::hash<int>{}(k.end2);
		return h1 ^ (h2 << 1U) ^ (h3 << 2U) ^ (h4 << 3U);
	}
};

struct Neighbor {
	string name;
	string orient;

	bool operator==(const Neighbor& other) const {
		return name == other.name && orient == other.orient;
	}
};

struct NeighborHash {
	size_t operator()(const Neighbor& n) const {
		size_t h1 = std::hash<string>{}(n.name);
		size_t h2 = std::hash<string>{}(n.orient);
		return h1 ^ (h2 << 1U);
	}
};

struct WrittenLink {
	string from;
	string from_orient;
	string to;
	string to_orient;

	bool operator==(const WrittenLink& other) const {
		return from == other.from && from_orient == other.from_orient && to == other.to && to_orient == other.to_orient;
	}
};

struct WrittenLinkHash {
	size_t operator()(const WrittenLink& w) const {
		size_t h1 = std::hash<string>{}(w.from);
		size_t h2 = std::hash<string>{}(w.from_orient);
		size_t h3 = std::hash<string>{}(w.to);
		size_t h4 = std::hash<string>{}(w.to_orient);
		return h1 ^ (h2 << 1U) ^ (h3 << 2U) ^ (h4 << 3U);
	}
};

vector<string> split_tab(const string& line) {
	vector<string> fields;
	std::stringstream ss(line);
	string field;
	while (std::getline(ss, field, '\t')) {
		fields.push_back(field);
	}
	return fields;
}

int parse_overlap1(const string& cigar, bool print_error_on_i) {
	int overlap1 = 0;
	string len_string;
	for (char c : cigar) {
		if (std::isdigit(static_cast<unsigned char>(c))) {
			len_string.push_back(c);
		} else {
			int n = len_string.empty() ? 0 : std::stoi(len_string);
			if (c == 'M') {
				overlap1 += n;
			} else if (c == 'D') {
				overlap1 += n;
			} else if (c == 'I') {
				if (print_error_on_i) {
					std::cout << "ERROR: bluntify only works with perfect overlaps" << std::endl;
				}
			}
			len_string.clear();
		}
	}
	return overlap1;
}

int parse_overlap2(const string& cigar, bool print_error_on_i) {
	int overlap2 = 0;
	string len_string;
	for (char c : cigar) {
		if (std::isdigit(static_cast<unsigned char>(c))) {
			len_string.push_back(c);
		} else {
			int n = len_string.empty() ? 0 : std::stoi(len_string);
			if (c == 'M') {
				overlap2 += n;
			} else if (c == 'I') {
				overlap2 += n;
				if (print_error_on_i) {
					std::cout << "ERROR: bluntify only works with perfect overlaps" << std::endl;
				}
			}
			len_string.clear();
		}
	}
	return overlap2;
}

string join_tail_with_tabs(const vector<string>& fields, int start) {
	if (start >= static_cast<int>(fields.size())) {
		return "";
	}
	string out;
	for (int i = start; i < static_cast<int>(fields.size()); ++i) {
		if (i > start) {
			out += "\t";
		}
		out += fields[i];
	}
	return out;
}

string join_tab(const vector<string>& fields) {
	return join_tail_with_tabs(fields, 0);
}

vector<string> split_name_parts(const string& s) {
	vector<string> parts;
	std::stringstream ss(s);
	string chunk;
	while (getline(ss, chunk, '_')) {
		parts.push_back(chunk);
	}
	return parts;
}

string orient_from_end1(int end1) {
	return string(1, "-+"[end1]);
}

string orient_from_end2(int end2) {
	return string(1, "+-"[end2]);
}

void parse_gfa_for_overlap_steps(
	const string& gfa_in,
	unordered_map<string, ContigSides>& list_of_contigs,
	vector<string>& contig_order,
	unordered_map<string, int>& length_of_contigs,
	unordered_map<string, std::streampos>& location_of_contigs_in_gfa,
	vector<Link>& list_of_links,
	unordered_map<string, string>* coverage,
	bool print_error_on_i
) {
	ifstream gfa_open(gfa_in);
	if (!gfa_open) {
		throw std::runtime_error("Cannot open input GFA: " + gfa_in);
	}

	unordered_set<LinkKey, LinkKeyHash> set_of_already_appended_links;

	while (true) {
		std::streampos pos = gfa_open.tellg();
		string line;
		if (!getline(gfa_open, line)) {
			break;
		}
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}

		if (line.rfind("S", 0) == 0) {
			vector<string> fields = split_tab(line);
			if (fields.size() < 2) {
				continue;
			}
			const string& name = fields[1];
			if (list_of_contigs.find(name) == list_of_contigs.end()) {
				list_of_contigs[name] = ContigSides{};
				contig_order.push_back(name);
			} else {
				list_of_contigs[name] = ContigSides{};
			}

			int length = 0;
			if (fields.size() > 2) {
				length = static_cast<int>(fields[2].size());
			}
			length_of_contigs[name] = length;
			location_of_contigs_in_gfa[name] = pos;

			if (coverage != nullptr) {
				for (const string& field : fields) {
					string upper = field;
					std::transform(upper.begin(), upper.end(), upper.begin(),
								   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
					if (upper.rfind("DP:", 0) == 0) {
						(*coverage)[name] = field;
						break;
					}
				}
			}
		} else if (line.rfind("L", 0) == 0) {
			vector<string> fields = split_tab(line);
			if (fields.size() < 6) {
				continue;
			}

			string name1 = fields[1];
			int end1 = static_cast<int>(fields[2] == "+");
			string name2 = fields[3];
			int end2 = static_cast<int>(fields[4] == "-");

			LinkKey key{name1, end1, name2, end2};
			if (set_of_already_appended_links.find(key) != set_of_already_appended_links.end()) {
				continue;
			}

			set_of_already_appended_links.insert(key);
			set_of_already_appended_links.insert(LinkKey{name2, end2, name1, end1});

			const string& cigar = fields[5];
			int overlap1 = parse_overlap1(cigar, print_error_on_i);
			int overlap2 = parse_overlap2(cigar, print_error_on_i);

			list_of_links.push_back(Link{name1, end1, overlap1, name2, end2, overlap2});
			int idx = static_cast<int>(list_of_links.size()) - 1;
			list_of_contigs[name1].left.reserve(list_of_contigs[name1].left.size());
			if (end1 == 0) {
				list_of_contigs[name1].left.push_back(idx);
			} else {
				list_of_contigs[name1].right.push_back(idx);
			}

			if (end2 == 0) {
				list_of_contigs[name2].left.push_back(idx);
			} else {
				list_of_contigs[name2].right.push_back(idx);
			}
		}
	}
}

}  // namespace

std::string reverse_complement(const std::string& seq) {
	unordered_map<char, char> complement{{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};
	string out;
	out.reserve(seq.size());
	for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
		auto found = complement.find(*it);
		out.push_back(found == complement.end() ? *it : found->second);
	}
	return out;
}

void basic_overlap_removal(const std::string& gfa_in, const std::string& gfa_out) {
	unordered_map<string, ContigSides> list_of_contigs;
	vector<string> contig_order;
	unordered_map<string, int> length_of_contigs;
	unordered_map<string, std::streampos> location_of_contigs_in_gfa;
	vector<Link> list_of_links;

	parse_gfa_for_overlap_steps(gfa_in, list_of_contigs, contig_order, length_of_contigs,
								location_of_contigs_in_gfa, list_of_links, nullptr, true);

	unordered_map<string, pair<int, int>> trimmed_lengths;
	for (const string& contig : contig_order) {
		int min_length_left = length_of_contigs[contig];
		int max_length_left = 0;
		if (list_of_contigs[contig].left.empty()) {
			min_length_left = 0;
		}

		int min_length_right = length_of_contigs[contig];
		int max_length_right = 0;
		if (list_of_contigs[contig].right.empty()) {
			min_length_right = 0;
		}

		for (int link_idx : list_of_contigs[contig].left) {
			int overlap_length = list_of_links[link_idx].overlap1;
			if (overlap_length > max_length_left) {
				max_length_left = overlap_length;
			}
			if (overlap_length < min_length_left) {
				min_length_left = overlap_length;
			}
		}

		for (int link_idx : list_of_contigs[contig].right) {
			int overlap_length = list_of_links[link_idx].overlap1;
			if (overlap_length > max_length_right) {
				max_length_right = overlap_length;
			}
			if (overlap_length < min_length_right) {
				min_length_right = overlap_length;
			}
		}

		int trim_left = std::min(min_length_left, length_of_contigs[contig] - max_length_right);
		int trim_right = std::min(min_length_right, length_of_contigs[contig] - max_length_left);

		for (int link_idx : list_of_contigs[contig].left) {
			list_of_links[link_idx].overlap1 -= trim_left;
			list_of_links[link_idx].overlap2 -= trim_left;
		}
		for (int link_idx : list_of_contigs[contig].right) {
			if (list_of_links[link_idx].name1 != list_of_links[link_idx].name2) {
				list_of_links[link_idx].overlap1 -= trim_right;
			}
		}

		trimmed_lengths[contig] = {trim_left, trim_right};
	}

	ofstream fo(gfa_out);
	ifstream fi(gfa_in);
	if (!fo || !fi) {
		throw std::runtime_error("Cannot open files for basic_overlap_removal output");
	}

	string line;
	while (getline(fi, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (!line.empty() && line[0] == 'S') {
			vector<string> ls = split_tab(line);
			if (ls.size() < 2) {
				continue;
			}
			string seq;
			if (ls.size() > 2) {
				seq = ls[2];
			}

			int left = trimmed_lengths[ls[1]].first;
			int right = length_of_contigs[ls[1]] - trimmed_lengths[ls[1]].second;
			if (left < 0) {
				left = 0;
			}
			if (right < left) {
				right = left;
			}
			string clipped = seq.substr(static_cast<size_t>(left), static_cast<size_t>(right - left));
			fo << "S\t" << ls[1] << "\t" << clipped << "\t" << join_tail_with_tabs(ls, 3) << "\n";
		}
	}

	for (const string& contig : contig_order) {
		for (int link_idx : list_of_contigs[contig].left) {
			if (list_of_links[link_idx].name1 != "None") {
				fo << "L\t" << list_of_links[link_idx].name1 << "\t"
				   << orient_from_end1(list_of_links[link_idx].end1) << "\t"
				   << list_of_links[link_idx].name2 << "\t"
				   << orient_from_end2(list_of_links[link_idx].end2) << "\t"
				   << list_of_links[link_idx].overlap1 << "M\n";
				list_of_links[link_idx].name1 = "None";
			}
		}

		for (int link_idx : list_of_contigs[contig].right) {
			if (list_of_links[link_idx].name1 != "None") {
				fo << "L\t" << list_of_links[link_idx].name1 << "\t"
				   << orient_from_end1(list_of_links[link_idx].end1) << "\t"
				   << list_of_links[link_idx].name2 << "\t"
				   << orient_from_end2(list_of_links[link_idx].end2) << "\t"
				   << list_of_links[link_idx].overlap1 << "M\n";
				list_of_links[link_idx].name1 = "None";
			}
		}
	}
}

void fancier_overlap_removal(const std::string& gfa_in, const std::string& gfa_out,
							  int short_contig_length) {
	unordered_map<string, ContigSides> list_of_contigs;
	vector<string> contig_order;
	unordered_map<string, int> length_of_contigs;
	unordered_map<string, std::streampos> location_of_contigs_in_gfa;
	vector<Link> list_of_links;
	unordered_map<string, string> coverage;

	parse_gfa_for_overlap_steps(gfa_in, list_of_contigs, contig_order, length_of_contigs,
								location_of_contigs_in_gfa, list_of_links, &coverage, false);

	unordered_map<string, ContigSides> new_list_of_contigs;
	vector<string> new_contig_order;
	vector<Link> new_list_of_links;
	unordered_set<string> contigs_already_dealt_with;
	unordered_map<string, pair<string, string>> old_contigs_to_new_contigs;

	auto set_contig_lists = [&](const string& name, const vector<int>& left, const vector<int>& right) {
		if (new_list_of_contigs.find(name) == new_list_of_contigs.end()) {
			new_contig_order.push_back(name);
		}
		new_list_of_contigs[name] = ContigSides{left, right};
	};

	for (const string& contig : contig_order) {
		int max_length_left = 0;
		int max_length_right = 0;
		vector<pair<int, int>> list_breakpoints_left;
		vector<pair<int, int>> list_breakpoints_right;

		for (int link_idx : list_of_contigs[contig].left) {
			int overlap_length = list_of_links[link_idx].overlap1;
			if (overlap_length > max_length_left) {
				max_length_left = overlap_length;
			}
			list_breakpoints_left.push_back({overlap_length, link_idx});
		}
		if (list_breakpoints_left.empty()) {
			list_breakpoints_left.push_back({0, -1});
		}

		for (int link_idx : list_of_contigs[contig].right) {
			int overlap_length = list_of_links[link_idx].overlap1;
			if (overlap_length > max_length_right) {
				max_length_right = overlap_length;
			}
			list_breakpoints_right.push_back({length_of_contigs[contig] - overlap_length, link_idx});
		}
		if (list_breakpoints_right.empty()) {
			list_breakpoints_right.push_back({length_of_contigs[contig], -1});
		}

		if (max_length_left + max_length_right < length_of_contigs[contig]) {
			contigs_already_dealt_with.insert(contig);

			vector<pair<int, int>> breakpoints = list_breakpoints_left;
			breakpoints.insert(breakpoints.end(), list_breakpoints_right.begin(), list_breakpoints_right.end());
			std::stable_sort(breakpoints.begin(), breakpoints.end(),
							 [](const pair<int, int>& a, const pair<int, int>& b) {
								 return a.first < b.first;
							 });

			string new_contig_name;
			for (int sub = 0; sub < static_cast<int>(breakpoints.size()); ++sub) {
				if (sub > 0 && breakpoints[sub - 1].first != breakpoints[sub].first) {
					new_contig_name = contig + "_" + std::to_string(breakpoints[sub - 1].first) + "_" +
									  std::to_string(breakpoints[sub].first);
					int n = 2;
					while (sub - n >= 0 && breakpoints[sub - 1].first == breakpoints[sub - n].first) {
						n += 1;
					}
					if (sub - n >= 0) {
						string past_contig_name = contig + "_" + std::to_string(breakpoints[sub - n].first) + "_" +
												  std::to_string(breakpoints[sub - 1].first);

						new_list_of_links.push_back(Link{new_contig_name, 0, 0, past_contig_name, 1, 0});
						set_contig_lists(new_contig_name,
										 vector<int>{static_cast<int>(new_list_of_links.size()) - 1},
										 vector<int>{});
					} else {
						set_contig_lists(new_contig_name, vector<int>{}, vector<int>{});
					}
				}

				if (sub == 0) {
					int n = 1;
					while (sub + n < static_cast<int>(breakpoints.size()) &&
						   std::to_string(breakpoints[sub].first) == std::to_string(breakpoints[sub + n].first)) {
						n += 1;
					}
					string future_contig_name = contig + "_" + std::to_string(breakpoints[sub].first) + "_" +
												std::to_string(breakpoints[sub + n].first);
					old_contigs_to_new_contigs[contig] = {future_contig_name, ""};
				}
				if (sub == static_cast<int>(breakpoints.size()) - 1) {
					old_contigs_to_new_contigs[contig].second = new_contig_name;
				}
			}

			for (int sub = 0; sub < static_cast<int>(breakpoints.size()); ++sub) {
				if (breakpoints[sub].second != -1) {
					Link link = list_of_links[breakpoints[sub].second];
					link.overlap1 = 0;
					link.overlap2 = 0;
					if (link.name1 == "-1") {
						continue;
					}

					if ((link.name1 == contig && link.end1 == 0) || (link.name2 == contig && link.end2 == 0)) {
						int n = 1;
						while (sub + n < static_cast<int>(breakpoints.size()) &&
							   std::to_string(breakpoints[sub].first) == std::to_string(breakpoints[sub + n].first)) {
							n += 1;
						}
						string future_contig_name = contig + "_" + std::to_string(breakpoints[sub].first) + "_" +
													std::to_string(breakpoints[sub + n].first);
						if (link.name1 == contig && link.end1 == 0) {
							link.name1 = future_contig_name;
						} else if (link.name2 == contig && link.end2 == 0) {
							link.name2 = future_contig_name;
						}
						new_list_of_links.push_back(link);
						new_list_of_contigs[future_contig_name].left.push_back(static_cast<int>(new_list_of_links.size()) - 1);

					} else {
						int n = 1;
						while (sub - n >= 0 && breakpoints[sub].first == breakpoints[sub - n].first) {
							n += 1;
						}
						string curr_contig_name = contig + "_" + std::to_string(breakpoints[sub - n].first) + "_" +
												  std::to_string(breakpoints[sub].first);
						if (link.name1 == contig && link.end1 == 1) {
							link.name1 = curr_contig_name;
						}
						if (link.name2 == contig && link.end2 == 1) {
							link.name2 = curr_contig_name;
						}
						new_list_of_links.push_back(link);
						new_list_of_contigs[curr_contig_name].right.push_back(static_cast<int>(new_list_of_links.size()) - 1);
					}

					list_of_links[breakpoints[sub].second].name1 = "-1";
					list_of_links[breakpoints[sub].second].name2 = "-1";
				}
			}

		} else {
			string new_name = contig + "_0_" + std::to_string(length_of_contigs[contig]);
			set_contig_lists(new_name, vector<int>{}, vector<int>{});
			old_contigs_to_new_contigs[contig] = {new_name, new_name};
		}
	}

	for (Link& link : new_list_of_links) {
		if (old_contigs_to_new_contigs.find(link.name1) != old_contigs_to_new_contigs.end()) {
			if (link.end1 == 0) {
				link.name1 = old_contigs_to_new_contigs[link.name1].first;
			} else {
				link.name1 = old_contigs_to_new_contigs[link.name1].second;
			}
		}
		if (old_contigs_to_new_contigs.find(link.name2) != old_contigs_to_new_contigs.end()) {
			if (link.end2 == 0) {
				link.name2 = old_contigs_to_new_contigs[link.name2].first;
			} else {
				link.name2 = old_contigs_to_new_contigs[link.name2].second;
			}
		}
	}

	ofstream fo(gfa_out);
	if (!fo) {
		throw std::runtime_error("Cannot open output GFA: " + gfa_out);
	}

	unordered_map<string, int> new_contigs_length;
	for (const string& contig : new_contig_order) {
		ifstream open2(gfa_in);
		if (!open2) {
			throw std::runtime_error("Cannot open input GFA: " + gfa_in);
		}

		vector<string> contig_parts = split_name_parts(contig);
		string original_name;
		for (int i = 0; i < static_cast<int>(contig_parts.size()) - 2; ++i) {
			if (i > 0) {
				original_name += "_";
			}
			original_name += contig_parts[i];
		}

		open2.seekg(location_of_contigs_in_gfa[original_name]);
		string sline;
		getline(open2, sline);
		if (!sline.empty() && sline.back() == '\r') {
			sline.pop_back();
		}
		vector<string> ls = split_tab(sline);
		string seq;
		if (ls.size() > 2) {
			seq = ls[2];
		}

		string left = contig_parts[contig_parts.size() - 2];
		string right = contig_parts[contig_parts.size() - 1];
		int ileft = std::stoi(left);
		int iright = std::stoi(right);

		if (iright - ileft >= short_contig_length) {
			fo << "S\t" << contig << "\t" << seq.substr(static_cast<size_t>(ileft), static_cast<size_t>(iright - ileft)) << "\t";
			if (coverage.find(original_name) != coverage.end()) {
				fo << coverage[original_name] << "\n";
			} else {
				fo << "\n";
			}
		}
		new_contigs_length[contig] = iright - ileft;
	}

	for (const string& contig : new_contig_order) {
		for (int link_idx : new_list_of_contigs[contig].left) {
			if (new_list_of_links[link_idx].name1 != "None") {
				if (new_list_of_links[link_idx].overlap1 == 0) {
					fo << "L\t" << new_list_of_links[link_idx].name1 << "\t"
					   << orient_from_end1(new_list_of_links[link_idx].end1) << "\t"
					   << new_list_of_links[link_idx].name2 << "\t"
					   << orient_from_end2(new_list_of_links[link_idx].end2) << "\t"
					   << new_list_of_links[link_idx].overlap1 << "M\n";
				}
				new_list_of_links[link_idx].name1 = "None";
			}
		}

		for (int link_idx : new_list_of_contigs[contig].right) {
			if (new_list_of_links[link_idx].name1 != "None") {
				if (new_list_of_links[link_idx].overlap1 == 0) {
					fo << "L\t" << new_list_of_links[link_idx].name1 << "\t"
					   << orient_from_end1(new_list_of_links[link_idx].end1) << "\t"
					   << new_list_of_links[link_idx].name2 << "\t"
					   << orient_from_end2(new_list_of_links[link_idx].end2) << "\t"
					   << new_list_of_links[link_idx].overlap1 << "M\n";
				}
				new_list_of_links[link_idx].name1 = "None";
			}
		}
	}
}

void remove_contigs_of_length_0(const std::string& gfa_in, const std::string& gfa_out) {
	unordered_map<string, string> contig_seqs;
	unordered_map<string, string> contig_lines;
	vector<string> contig_order;
	vector<vector<string>> links;

	{
		ifstream fi(gfa_in);
		if (!fi) {
			throw std::runtime_error("Cannot open input GFA: " + gfa_in);
		}

		string line;
		while (getline(fi, line)) {
			if (!line.empty() && line.back() == '\r') {
				line.pop_back();
			}
			if (line.rfind("S", 0) == 0) {
				vector<string> fields = split_tab(line);
				if (fields.size() < 2) {
					continue;
				}
				string name = fields[1];
				string seq;
				if (fields.size() > 2) {
					seq = fields[2];
				}
				if (contig_seqs.find(name) == contig_seqs.end()) {
					contig_order.push_back(name);
				}
				contig_seqs[name] = seq;
				contig_lines[name] = line + "\n";
			} else if (line.rfind("L", 0) == 0) {
				links.push_back(split_tab(line));
			}
		}
	}

	unordered_set<string> zero_length_contigs;
	for (const auto& kv : contig_seqs) {
		if (kv.second.empty()) {
			zero_length_contigs.insert(kv.first);
		}
	}

	unordered_map<string, unordered_set<Neighbor, NeighborHash>> left_neighbors;
	unordered_map<string, unordered_set<Neighbor, NeighborHash>> right_neighbors;
	for (const auto& kv : contig_seqs) {
		left_neighbors[kv.first] = {};
		right_neighbors[kv.first] = {};
	}

	for (const vector<string>& fields : links) {
		if (fields.size() < 5) {
			continue;
		}
		string from_name = fields[1];
		string to_name = fields[3];
		string from_orient = fields[2];
		string to_orient = fields[4];
		if (from_orient == "+") {
			right_neighbors[from_name].insert(Neighbor{to_name, (to_orient == "+") ? "-" : "+"});
		} else {
			left_neighbors[from_name].insert(Neighbor{to_name, (to_orient == "+") ? "-" : "+"});
		}
		if (to_orient == "+") {
			left_neighbors[to_name].insert(Neighbor{from_name, from_orient});
		} else {
			right_neighbors[to_name].insert(Neighbor{from_name, from_orient});
		}
	}

	vector<string> zero_contigs_order(zero_length_contigs.begin(), zero_length_contigs.end());
	for (const string& zero_contig : zero_contigs_order) {
		vector<Neighbor> lefts(left_neighbors[zero_contig].begin(), left_neighbors[zero_contig].end());
		vector<Neighbor> rights(right_neighbors[zero_contig].begin(), right_neighbors[zero_contig].end());

		for (const Neighbor& left_neighbor : lefts) {
			for (const Neighbor& right_neighbor : rights) {
				if (left_neighbor.name == zero_contig || right_neighbor.name == zero_contig) {
					continue;
				}

				if (left_neighbor.orient == "+") {
					right_neighbors[left_neighbor.name].insert(Neighbor{right_neighbor.name, right_neighbor.orient});
				} else {
					left_neighbors[left_neighbor.name].insert(Neighbor{right_neighbor.name, right_neighbor.orient});
				}

				if (left_neighbor.name != right_neighbor.name || left_neighbor.orient != right_neighbor.orient) {
					if (right_neighbor.orient == "-") {
						left_neighbors[right_neighbor.name].insert(Neighbor{left_neighbor.name, left_neighbor.orient});
					} else {
						right_neighbors[right_neighbor.name].insert(Neighbor{left_neighbor.name, left_neighbor.orient});
					}
				}
			}
		}
	}

	{
		ofstream fo(gfa_out);
		if (!fo) {
			throw std::runtime_error("Cannot open output GFA: " + gfa_out);
		}
		for (const string& name : contig_order) {
			if (zero_length_contigs.find(name) == zero_length_contigs.end()) {
				fo << contig_lines[name];
			}
		}
	}

	{
		ofstream fo(gfa_out, std::ios::app);
		if (!fo) {
			throw std::runtime_error("Cannot append to output GFA: " + gfa_out);
		}
		unordered_set<WrittenLink, WrittenLinkHash> written_links;

		for (const string& from_name : contig_order) {
			if (zero_length_contigs.find(from_name) != zero_length_contigs.end()) {
				continue;
			}

			for (const Neighbor& to : right_neighbors[from_name]) {
				if (zero_length_contigs.find(to.name) != zero_length_contigs.end()) {
					continue;
				}
				WrittenLink link_tuple{from_name, "+", to.name, to.orient};
				if (written_links.find(link_tuple) == written_links.end()) {
					fo << "L\t" << from_name << "\t+\t" << to.name << "\t"
					   << ((to.orient == "+") ? "-" : "+") << "\t0M\n";
					written_links.insert(link_tuple);
					written_links.insert(WrittenLink{to.name, to.orient, from_name, "+"});
				}
			}

			for (const Neighbor& to : left_neighbors[from_name]) {
				if (zero_length_contigs.find(to.name) != zero_length_contigs.end()) {
					continue;
				}
				WrittenLink link_tuple{from_name, "-", to.name, to.orient};
				if (written_links.find(link_tuple) == written_links.end()) {
					fo << "L\t" << from_name << "\t-\t" << to.name << "\t"
					   << ((to.orient == "+") ? "-" : "+") << "\t0M\n";
					written_links.insert(link_tuple);
					written_links.insert(WrittenLink{to.name, to.orient, from_name, "-"});
				}
			}
		}
	}
}

/**
 * @brief Processes a GFA file to remove overlaps and isolated contigs, producing a "bluntified" output.
 *
 * @param input Path to the input GFA file.
 * @param output Path where the processed (bluntified) GFA file will be written.
 * @param int trim_isolated_contigs_length If positive, ends of isolated contigs will be trimmed (the idea being that these kmers are already elsewhere in non-isolated contigs).
 * @param tmpdir Directory to use for storing temporary files during processing.
 */
void bluntify(const std::string& input, const std::string& output,
			  int trim_isolated_contigs_length, const std::string& tmpdir) {

	string intermediate_gfa = tmpdir + "/intermediate_gfa.tmp.gfa";
    string intermediate_gfa_2 = tmpdir + "/intermediate_gfa_2.tmp.gfa";
	basic_overlap_removal(input, intermediate_gfa);
	fancier_overlap_removal(intermediate_gfa, intermediate_gfa + ".fancy.gfa");
	remove_contigs_of_length_0(intermediate_gfa + ".fancy.gfa", intermediate_gfa_2);

	if (trim_isolated_contigs_length > 0){
        //go through all contigs and trim ends of isolated contigs
        unordered_set<string> isolated_contigs;
        {
            ifstream fi(intermediate_gfa_2);

            string line;
            while (getline(fi, line)) {
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                if (line.rfind("S", 0) == 0) {
                    vector<string> fields = split_tab(line);
                    if (fields.size() < 2) {
                        continue;
                    }
                    string name = fields[1];
                    isolated_contigs.insert(name);
                } else if (line.rfind("L", 0) == 0) {
                    vector<string> fields = split_tab(line);
                    if (fields.size() < 4) {
                        continue;
                    }
                    string from_name = fields[1];
                    string to_name = fields[3];
                    isolated_contigs.erase(from_name);
                    isolated_contigs.erase(to_name);
                }
            }
        }
        //write output with isolated contigs trimmed of trim_isolated_contigs_length from each end
        ofstream fo(output);
        if (!fo) {
            throw std::runtime_error("Cannot open output GFA: " + output);
        }
        ifstream fi(intermediate_gfa_2);
        string line;
        while (getline(fi, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (line.rfind("S", 0) == 0) {
                vector<string> fields = split_tab(line);
                if (fields.size() < 2) {
                    continue;
                }
                string name = fields[1];
                string seq;
                if (fields.size() > 2) {
                    seq = fields[2];
                }
                if (isolated_contigs.find(name) != isolated_contigs.end()) {
                    int trim_length = std::min(trim_isolated_contigs_length, static_cast<int>(seq.size()) / 2);
                    seq = seq.substr(static_cast<size_t>(trim_length), static_cast<size_t>(seq.size() - 2 * trim_length));
                    fields[2] = seq;
                    fo << join_tab(fields) << "\n";
                } else {
                    fo << line << "\n";
                }
            } else {
                fo << line << "\n";
            }
        }
	}
	else {
		std::filesystem::copy_file(
			intermediate_gfa_2,
			output,
			std::filesystem::copy_options::overwrite_existing
		);
	}

	std::filesystem::remove(intermediate_gfa);
	std::filesystem::remove(intermediate_gfa_2);
	std::filesystem::remove(intermediate_gfa + ".fancy.gfa");
}

int bluntify_main(int argc, char** argv) {
	if (argc == 2 && string(argv[1]) == "--version") {
		std::cout << argv[0] << " " << kVersion << "\n";
		return 0;
	}

	bool no_overlaps = false;
	string tmpdir = ".";
	vector<string> positional;

	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];
		if (arg == "-n" || arg == "--no_overlaps") {
			no_overlaps = true;
		} else if (arg == "--tmpdir") {
			if (i + 1 >= argc) {
				std::cerr << "Error: --tmpdir requires a value\n";
				return 1;
			}
			tmpdir = argv[++i];
		} else if (!arg.empty() && arg[0] == '-') {
			std::cerr << "Error: unknown option " << arg << "\n";
			return 1;
		} else {
			positional.push_back(arg);
		}
	}

	if (positional.size() != 2) {
		std::cerr << "Usage: " << argv[0]
				  << " input output [-n|--no_overlaps] [--tmpdir TMPDIR] [--version]\n";
		return 1;
	}

	bluntify(positional[0], positional[1], no_overlaps, tmpdir);
	return 0;
}
