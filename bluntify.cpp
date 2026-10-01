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
 * - split_tab, join_tail_with_tabs: String manipulation helpers for parsing GFA lines.
 * - parse_overlaps: Parse CIGAR strings to determine overlap lengths.
 * - orient_from_end1, orient_from_end2: Determine orientation symbols from link ends.
 * - parse_gfa_for_overlap_steps: Parses a GFA file and populates data structures for contigs and links.
 *
 * Self-links are supported: a circular self-link (A + A + ov) is trimmed once (the circle keeps len-ov bases),
 * a hairpin (A + A - ov, the end of A is a palindrome of length ov) is trimmed by ov/2.
 *
 * Usage:
 *   bluntify <input.gfa> <output.gfa> [-t|--trim_isolated LENGTH] [--tmpdir TMPDIR] [--version]
 *
 * Dependencies:
 *   Requires C++17 or later for <filesystem> and other standard library features.
 *
 * Author: Roland Faure
 * Version: 0.2
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

const string kVersion = "0.2";

struct Link {
	string name1;
	int end1 = 0; //0 = left end of name1, 1 = right end
	int overlap1 = 0;
	string name2;
	int end2 = 0; //0 = left end of name2, 1 = right end
	int overlap2 = 0;
	bool written = false; //the link has already been output
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

/**
 * @brief Reverse complement, kept internal to bluntify (and named differently) so that it can never be confused with
 * reverse_complement(string&) of basic_graph_manipulation.h.
 */
[[maybe_unused]] string bluntify_reverse_complement(const string& seq) {
	string out;
	out.reserve(seq.size());
	for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
		switch (*it) {
			case 'A': out.push_back('T'); break;
			case 'C': out.push_back('G'); break;
			case 'G': out.push_back('C'); break;
			case 'T': out.push_back('A'); break;
			default: out.push_back(*it); break;
		}
	}
	return out;
}

//same result as reading the fields with getline(..., '\t') (a trailing empty field is dropped), without a stringstream
vector<string> split_tab(const string& line) {
	vector<string> fields;
	size_t start = 0;
	while (start < line.size()) {
		size_t tab = line.find('\t', start);
		if (tab == string::npos) {
			fields.push_back(line.substr(start));
			break;
		}
		fields.push_back(line.substr(start, tab - start));
		start = tab + 1;
	}
	return fields;
}

/**
 * @brief Parses a CIGAR string. overlap1 is the length of the overlap on the first contig, overlap2 on the second.
 * M, = and X count on both contigs, D only on the first, I only on the second. "*" means no overlap (0).
 *
 * @return true if the overlap is not a perfect match (contains I, D or another operation)
 */
bool parse_overlaps(const string& cigar, int& overlap1, int& overlap2) {
	overlap1 = 0;
	overlap2 = 0;
	bool imperfect = false;
	if (cigar == "*") {
		return false;
	}
	string len_string;
	for (char c : cigar) {
		if (std::isdigit(static_cast<unsigned char>(c))) {
			len_string.push_back(c);
		} else {
			int n = len_string.empty() ? 0 : std::stoi(len_string);
			if (c == 'M' || c == '=' || c == 'X') {
				overlap1 += n;
				overlap2 += n;
			} else if (c == 'D') {
				overlap1 += n;
				imperfect = true;
			} else if (c == 'I') {
				overlap2 += n;
				imperfect = true;
			} else {
				imperfect = true;
			}
			len_string.clear();
		}
	}
	return imperfect;
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

string orient_from_end1(int end1) {
	return string(1, "-+"[end1]);
}

string orient_from_end2(int end2) {
	return string(1, "+-"[end2]);
}

//circular self-link, e.g. A + A +: it is both on the left and on the right of the contig
bool is_circular_self_link(const Link& link) {
	return link.name1 == link.name2 && link.end1 != link.end2;
}

//hairpin, e.g. A + A -: the end of the contig is a palindrome of length overlap, it appears twice in the same side
bool is_hairpin(const Link& link) {
	return link.name1 == link.name2 && link.end1 == link.end2;
}

/**
 * @brief Number of bases that the contig must lose at this end to make the link blunt: the overlap, except for hairpins where
 * trimming t bases at the end of the contig decreases the overlap by 2t
 */
int effective_overlap(const Link& link) {
	return is_hairpin(link) ? link.overlap1 / 2 : link.overlap1;
}

void write_link(ofstream& fo, const Link& link) {
	fo << "L\t" << link.name1 << "\t" << orient_from_end1(link.end1) << "\t"
	   << link.name2 << "\t" << orient_from_end2(link.end2) << "\t"
	   << link.overlap1 << "M\n";
}

/**
 * @brief Parses the GFA. Contigs are listed in contig_order in the order of their S lines (also when L lines come first);
 * links pointing to contigs without S line are reported in contigs_without_S_line.
 */
void parse_gfa_for_overlap_steps(
	const string& gfa_in,
	unordered_map<string, ContigSides>& list_of_contigs,
	vector<string>& contig_order,
	unordered_map<string, int>& length_of_contigs,
	unordered_map<string, std::streampos>& location_of_contigs_in_gfa,
	vector<Link>& list_of_links,
	unordered_map<string, string>* coverage,
	bool print_error_on_imperfect_overlaps
) {
	ifstream gfa_open(gfa_in);
	if (!gfa_open) {
		throw std::runtime_error("Cannot open input GFA: " + gfa_in);
	}

	unordered_set<LinkKey, LinkKeyHash> set_of_already_appended_links;
	unordered_set<string> contigs_with_S_line;
	int number_of_imperfect_overlaps = 0;

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
			//do not reset the entry: links of this contig may have been parsed before its S line
			list_of_contigs[name];
			if (contigs_with_S_line.insert(name).second) {
				contig_order.push_back(name);
			}

			int length = 0;
			if (fields.size() > 2) {
				length = static_cast<int>(fields[2].size());
			}
			length_of_contigs[name] = length;
			location_of_contigs_in_gfa[name] = pos;

			if (coverage != nullptr) {
				for (const string& field : fields) {
					if (field.size() >= 3 && std::toupper(static_cast<unsigned char>(field[0])) == 'D' &&
						std::toupper(static_cast<unsigned char>(field[1])) == 'P' && field[2] == ':') {
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

			const string& name1 = fields[1];
			int end1 = static_cast<int>(fields[2] == "+");
			const string& name2 = fields[3];
			int end2 = static_cast<int>(fields[4] == "-");

			LinkKey key{name1, end1, name2, end2};
			if (set_of_already_appended_links.find(key) != set_of_already_appended_links.end()) {
				continue;
			}

			set_of_already_appended_links.insert(key);
			set_of_already_appended_links.insert(LinkKey{name2, end2, name1, end1});

			int overlap1 = 0;
			int overlap2 = 0;
			if (parse_overlaps(fields[5], overlap1, overlap2)) {
				number_of_imperfect_overlaps += 1;
			}

			list_of_links.push_back(Link{name1, end1, overlap1, name2, end2, overlap2});
			int idx = static_cast<int>(list_of_links.size()) - 1;
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

	if (print_error_on_imperfect_overlaps && number_of_imperfect_overlaps > 0) {
		std::cout << "ERROR: bluntify only works with perfect overlaps, " << number_of_imperfect_overlaps
				  << " links of " << gfa_in << " have a CIGAR that is not only made of M/=/X" << std::endl;
	}
	if (print_error_on_imperfect_overlaps && contigs_with_S_line.size() != list_of_contigs.size()) {
		std::cout << "WARNING: in bluntify, " << list_of_contigs.size() - contigs_with_S_line.size()
				  << " contigs of " << gfa_in << " appear in L lines but have no S line, their links are ignored" << std::endl;
	}
}

}  // namespace

void basic_overlap_removal(const std::string& gfa_in, const std::string& gfa_out) {
	unordered_map<string, ContigSides> list_of_contigs;
	vector<string> contig_order;
	unordered_map<string, int> length_of_contigs;
	unordered_map<string, std::streampos> location_of_contigs_in_gfa;
	vector<Link> list_of_links;

	parse_gfa_for_overlap_steps(gfa_in, list_of_contigs, contig_order, length_of_contigs,
								location_of_contigs_in_gfa, list_of_links, nullptr, true);

	unordered_map<string, pair<int, int>> trimmed_lengths;
	int number_of_odd_hairpins = 0;
	for (const string& contig : contig_order) {
		const ContigSides& sides = list_of_contigs[contig];
		int length = length_of_contigs[contig];

		//trim the left end by the smallest overlap on the left, without eating into the overlaps on the right
		//(circular self-links are not counted on the right: trimming the left end shortens them too)
		int min_length_left = sides.left.empty() ? 0 : length;
		int max_length_right = 0;
		for (int link_idx : sides.left) {
			min_length_left = std::min(min_length_left, effective_overlap(list_of_links[link_idx]));
		}
		for (int link_idx : sides.right) {
			if (!is_circular_self_link(list_of_links[link_idx])) {
				max_length_right = std::max(max_length_right, list_of_links[link_idx].overlap1);
			}
		}
		int trim_left = std::max(0, std::min(min_length_left, length - max_length_right));
		for (int link_idx : sides.left) {
			list_of_links[link_idx].overlap1 -= trim_left;
			list_of_links[link_idx].overlap2 -= trim_left;
		}

		//then the right end, recomputing the overlaps after the left trimming (a circular self-link has just been shortened
		//by trim_left and must not be trimmed a second time)
		int min_length_right = sides.right.empty() ? 0 : length - trim_left;
		int max_length_left = 0;
		for (int link_idx : sides.right) {
			min_length_right = std::min(min_length_right, effective_overlap(list_of_links[link_idx]));
		}
		for (int link_idx : sides.left) {
			if (!is_circular_self_link(list_of_links[link_idx])) {
				max_length_left = std::max(max_length_left, list_of_links[link_idx].overlap1);
			}
		}
		int trim_right = std::max(0, std::min(min_length_right, length - trim_left - max_length_left));
		for (int link_idx : sides.right) {
			list_of_links[link_idx].overlap1 -= trim_right;
			list_of_links[link_idx].overlap2 -= trim_right;
		}

		for (int link_idx : sides.left) {
			if (is_hairpin(list_of_links[link_idx]) && list_of_links[link_idx].overlap1 % 2 == 1) {
				number_of_odd_hairpins += 1;
			}
		}
		for (int link_idx : sides.right) {
			if (is_hairpin(list_of_links[link_idx]) && list_of_links[link_idx].overlap1 % 2 == 1) {
				number_of_odd_hairpins += 1;
			}
		}

		trimmed_lengths[contig] = {trim_left, trim_right};
	}
	if (number_of_odd_hairpins > 0) {
		std::cout << "WARNING: in bluntify, " << number_of_odd_hairpins / 2
				  << " hairpin links have an odd overlap (not a palindrome), they cannot be made perfectly blunt" << std::endl;
	}

	ofstream fo(gfa_out);
	ifstream fi(gfa_in);
	if (!fo || !fi) {
		throw std::runtime_error("Cannot open files for basic_overlap_removal output");
	}

	unordered_set<string> contigs_already_written;
	string line;
	while (getline(fi, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (!line.empty() && line[0] == 'S') {
			vector<string> ls = split_tab(line);
			if (ls.size() < 2 || !contigs_already_written.insert(ls[1]).second) {
				continue;
			}
			const pair<int, int>& trimmed = trimmed_lengths[ls[1]];
			int length = length_of_contigs[ls[1]];
			int left = std::max(0, trimmed.first);
			int right = std::max(left, length - trimmed.second);
			fo << "S\t" << ls[1] << "\t";
			if (ls.size() > 2) {
				fo.write(ls[2].data() + left, right - left);
			}
			fo << "\t" << join_tail_with_tabs(ls, 3) << "\n";
		}
	}

	for (const string& contig : contig_order) {
		for (int link_idx : list_of_contigs[contig].left) {
			if (!list_of_links[link_idx].written) {
				write_link(fo, list_of_links[link_idx]);
				list_of_links[link_idx].written = true;
			}
		}

		for (int link_idx : list_of_contigs[contig].right) {
			if (!list_of_links[link_idx].written) {
				write_link(fo, list_of_links[link_idx]);
				list_of_links[link_idx].written = true;
			}
		}
	}
}

/**
 * @brief Splits the contigs at the breakpoints of the remaining overlaps and rewires the links as 0M.
 *
 * The overlap of each link is "cut" between its two sides: cut1 bases are removed on the side of name1 and cut2 on the side
 * of name2, with cut1 + cut2 = overlap. A side with a cut c is attached to the sub-contig starting at c (left end) or ending
 * at length-c (right end), so the overlapping sequence is kept exactly once on every walk. A contig can take the cuts only if
 * (max cut on the left) + (max cut on the right) <= length (when equal, a zero-length sub-contig is put at the junction and
 * removed later by remove_contigs_of_length_0).
 * - Contigs where this holds for all their overlaps ("splittable") take the whole overlap of the links they are the first
 *   (in GFA order) to see, as before.
 * - Links between two non-splittable contigs are cut on a non-splittable contig that can take all of them, if any (iterated),
 *   then by splitting the overlap between the two sides according to the budget of each contig end.
 * - Links that still cannot be made blunt are output with their remaining overlap and a warning is printed.
 */
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

	//rank of each contig in the GFA and whether all its overlaps can be cut on it
	unordered_map<string, int> rank_of_contig;
	unordered_map<string, bool> splittable;
	for (int r = 0; r < static_cast<int>(contig_order.size()); ++r) {
		const string& contig = contig_order[r];
		rank_of_contig[contig] = r;
		int max_length_left = 0;
		int max_length_right = 0;
		for (int link_idx : list_of_contigs[contig].left) {
			max_length_left = std::max(max_length_left, effective_overlap(list_of_links[link_idx]));
		}
		for (int link_idx : list_of_contigs[contig].right) {
			max_length_right = std::max(max_length_right, effective_overlap(list_of_links[link_idx]));
		}
		splittable[contig] = max_length_left + max_length_right <= length_of_contigs[contig];
	}

	//decide the cuts. -1 = not decided (yet)
	const int number_of_links = static_cast<int>(list_of_links.size());
	vector<int> cut1(number_of_links, -1);
	vector<int> cut2(number_of_links, -1);
	vector<bool> link_is_valid(number_of_links, true); //both contigs have an S line
	vector<int> pending; //links between two non-splittable contigs
	int number_of_links_to_missing_contigs = 0;
	for (int link_idx = 0; link_idx < number_of_links; ++link_idx) {
		const Link& link = list_of_links[link_idx];
		if (rank_of_contig.find(link.name1) == rank_of_contig.end() || rank_of_contig.find(link.name2) == rank_of_contig.end()) {
			link_is_valid[link_idx] = false;
			number_of_links_to_missing_contigs += 1;
			continue;
		}
		int overlap = effective_overlap(link);
		if (overlap <= 0) {
			cut1[link_idx] = 0;
			cut2[link_idx] = 0;
		} else if (is_hairpin(link)) {
			//both sides are the same end: it loses overlap/2 bases
			if (splittable[link.name1]) {
				cut1[link_idx] = overlap;
				cut2[link_idx] = overlap;
			} else {
				pending.push_back(link_idx);
			}
		} else {
			bool first_is_name1 = rank_of_contig[link.name1] <= rank_of_contig[link.name2];
			const string& first = first_is_name1 ? link.name1 : link.name2;
			const string& second = first_is_name1 ? link.name2 : link.name1;
			if (splittable[first]) {
				(first_is_name1 ? cut1 : cut2)[link_idx] = overlap;
				(first_is_name1 ? cut2 : cut1)[link_idx] = 0;
			} else if (splittable[second]) {
				(first_is_name1 ? cut2 : cut1)[link_idx] = overlap;
				(first_is_name1 ? cut1 : cut2)[link_idx] = 0;
			} else {
				pending.push_back(link_idx);
			}
		}
	}
	if (number_of_links_to_missing_contigs > 0) {
		std::cout << "WARNING: in bluntify, " << number_of_links_to_missing_contigs
				  << " links point to contigs without S line, they are dropped" << std::endl;
	}

	//links between non-splittable contigs. First, a non-splittable contig that can take all its pending overlaps takes them
	//(which frees its neighbours), until no such contig is left
	if (!pending.empty()) {
		unordered_map<string, vector<int>> pending_links_of_contig;
		for (int link_idx : pending) {
			pending_links_of_contig[list_of_links[link_idx].name1].push_back(link_idx);
			if (list_of_links[link_idx].name2 != list_of_links[link_idx].name1) {
				pending_links_of_contig[list_of_links[link_idx].name2].push_back(link_idx);
			}
		}
		//max pending overlap on the left / right of the contig
		auto pending_needs = [&](const string& contig) {
			pair<int, int> needs = {0, 0};
			for (int link_idx : pending_links_of_contig[contig]) {
				if (cut1[link_idx] != -1) {
					continue;
				}
				const Link& link = list_of_links[link_idx];
				int overlap = effective_overlap(link);
				if (link.name1 == contig) {
					(link.end1 == 0 ? needs.first : needs.second) = std::max(link.end1 == 0 ? needs.first : needs.second, overlap);
				}
				if (link.name2 == contig) {
					(link.end2 == 0 ? needs.first : needs.second) = std::max(link.end2 == 0 ? needs.first : needs.second, overlap);
				}
			}
			return needs;
		};

		vector<string> queue;
		for (const auto& kv : pending_links_of_contig) {
			queue.push_back(kv.first);
		}
		std::sort(queue.begin(), queue.end(), [&rank_of_contig](const string& a, const string& b) {
			return rank_of_contig[a] < rank_of_contig[b];
		});
		unordered_set<string> in_queue(queue.begin(), queue.end());
		while (!queue.empty()) {
			string contig = queue.back();
			queue.pop_back();
			in_queue.erase(contig);
			pair<int, int> needs = pending_needs(contig);
			if (needs.first + needs.second > length_of_contigs[contig]) {
				continue;
			}
			for (int link_idx : pending_links_of_contig[contig]) {
				if (cut1[link_idx] != -1) {
					continue;
				}
				const Link& link = list_of_links[link_idx];
				int overlap = effective_overlap(link);
				if (link.name1 == link.name2) { //self-link: cut on the left end for a circle, on its only end for a hairpin
					cut1[link_idx] = (is_hairpin(link) || link.end1 == 0) ? overlap : 0;
					cut2[link_idx] = (is_hairpin(link) || link.end2 == 0) ? overlap : 0;
				} else {
					cut1[link_idx] = (link.name1 == contig) ? overlap : 0;
					cut2[link_idx] = (link.name2 == contig) ? overlap : 0;
					const string& neighbor = (link.name1 == contig) ? link.name2 : link.name1;
					if (in_queue.insert(neighbor).second) {
						queue.push_back(neighbor);
					}
				}
			}
		}

		//then share the remaining overlaps between both sides, each contig end having a budget proportional to its needs
		unordered_map<string, pair<int, int>> budgets;
		for (const auto& kv : pending_links_of_contig) {
			pair<int, int> needs = pending_needs(kv.first);
			int length = length_of_contigs[kv.first];
			if (needs.first + needs.second <= length) {
				budgets[kv.first] = needs;
			} else {
				int left = static_cast<int>(static_cast<long long>(length) * needs.first / (needs.first + needs.second));
				//hairpins can only be cut on this contig: keep their part of the budget if possible
				pair<int, int> hairpin_needs = {0, 0};
				for (int link_idx : kv.second) {
					const Link& link = list_of_links[link_idx];
					if (cut1[link_idx] == -1 && is_hairpin(link)) {
						int& need = (link.end1 == 0) ? hairpin_needs.first : hairpin_needs.second;
						need = std::max(need, effective_overlap(link));
					}
				}
				if (hairpin_needs.first + hairpin_needs.second <= length) {
					left = std::min(std::max(left, hairpin_needs.first), length - hairpin_needs.second);
				}
				budgets[kv.first] = {left, length - left};
			}
		}
		for (int link_idx : pending) {
			if (cut1[link_idx] != -1) {
				continue;
			}
			const Link& link = list_of_links[link_idx];
			int overlap = effective_overlap(link);
			const pair<int, int>& budget1 = budgets[link.name1];
			const pair<int, int>& budget2 = budgets[link.name2];
			int available1 = (link.end1 == 0) ? budget1.first : budget1.second;
			int available2 = (link.end2 == 0) ? budget2.first : budget2.second;
			if (is_hairpin(link)) {
				if (overlap <= available1) {
					cut1[link_idx] = overlap;
					cut2[link_idx] = overlap;
				}
			} else if (overlap <= available1 + available2) {
				cut1[link_idx] = std::min(overlap, available1);
				cut2[link_idx] = overlap - cut1[link_idx];
			}
		}

		//last chance for the links left: use what the contigs have not used of their length (hairpins first, they have only one side)
		unordered_map<string, pair<int, int>> used; //max cut on the left / right of each contig
		auto use = [&used](const string& contig, int end, int cut) {
			int& u = (end == 0) ? used[contig].first : used[contig].second;
			u = std::max(u, cut);
		};
		vector<int> links_left;
		for (int link_idx : pending) {
			const Link& link = list_of_links[link_idx];
			if (cut1[link_idx] != -1) {
				use(link.name1, link.end1, cut1[link_idx]);
				use(link.name2, link.end2, cut2[link_idx]);
			} else {
				links_left.push_back(link_idx);
			}
		}
		std::stable_partition(links_left.begin(), links_left.end(),
							  [&list_of_links](int link_idx) { return is_hairpin(list_of_links[link_idx]); });
		for (int link_idx : links_left) {
			const Link& link = list_of_links[link_idx];
			int overlap = effective_overlap(link);
			pair<int, int> used1 = used[link.name1];
			int available1 = length_of_contigs[link.name1] - ((link.end1 == 0) ? used1.second : used1.first);
			int c1 = is_hairpin(link) ? overlap : std::min(overlap, std::max(0, available1));
			int c2 = is_hairpin(link) ? overlap : overlap - c1;
			//check the constraint of both contigs with the new cuts
			pair<int, int> new_used1 = used[link.name1];
			(link.end1 == 0 ? new_used1.first : new_used1.second) = std::max(link.end1 == 0 ? new_used1.first : new_used1.second, c1);
			pair<int, int> new_used2 = (link.name2 == link.name1) ? new_used1 : used[link.name2];
			(link.end2 == 0 ? new_used2.first : new_used2.second) = std::max(link.end2 == 0 ? new_used2.first : new_used2.second, c2);
			if (link.name2 == link.name1) {
				new_used1 = new_used2;
			}
			if (new_used1.first + new_used1.second <= length_of_contigs[link.name1] &&
				new_used2.first + new_used2.second <= length_of_contigs[link.name2]) {
				cut1[link_idx] = c1;
				cut2[link_idx] = c2;
				used[link.name1] = new_used1;
				used[link.name2] = new_used2;
			}
		}
	}

	struct Attachment {
		int position;
		int link_idx;
		int end; //end of the contig on which the link is
		int side; //1 if the link is attached by name1, 2 by name2
	};
	struct SubContig {
		string original_name;
		int start;
		int end;
	};

	unordered_map<string, ContigSides> new_list_of_contigs;
	vector<string> new_contig_order;
	vector<SubContig> new_contig_coordinates; //parallel to new_contig_order
	vector<Link> new_list_of_links;
	vector<int> new_index_of_link(number_of_links, -1);
	vector<int> number_of_sides_rewired(number_of_links, 0);

	auto sub_contig_name = [](const string& contig, int start, int end) {
		return contig + "_" + std::to_string(start) + "_" + std::to_string(end);
	};

	for (const string& contig : contig_order) {
		int length = length_of_contigs[contig];

		//where each link is attached on the contig
		vector<Attachment> attachments;
		for (int end = 0; end < 2; ++end) {
			const vector<int>& links_of_this_end = (end == 0) ? list_of_contigs[contig].left : list_of_contigs[contig].right;
			for (int i = 0; i < static_cast<int>(links_of_this_end.size()); ++i) {
				int link_idx = links_of_this_end[i];
				if (!link_is_valid[link_idx] || (i > 0 && links_of_this_end[i - 1] == link_idx)) { //hairpins are listed twice
					continue;
				}
				const Link& link = list_of_links[link_idx];
				for (int side = 1; side <= 2; ++side) {
					const string& name = (side == 1) ? link.name1 : link.name2;
					int link_end = (side == 1) ? link.end1 : link.end2;
					if (name != contig || link_end != end) {
						continue;
					}
					int cut = std::max(0, (side == 1) ? cut1[link_idx] : cut2[link_idx]); //undecided links are attached to the extremities
					attachments.push_back({(end == 0) ? cut : length - cut, link_idx, end, side});
				}
			}
		}
		std::stable_sort(attachments.begin(), attachments.end(),
						 [](const Attachment& a, const Attachment& b) {
							 return a.position < b.position;
						 });

		//the breakpoints delimit the sub-contigs. A breakpoint where links are attached both on the left and on the right
		//gets a zero-length sub-contig, so that links entering there can exit there
		unordered_set<int> left_positions;
		unordered_set<int> right_positions;
		vector<int> positions = {0, length};
		for (const Attachment& attachment : attachments) {
			positions.push_back(attachment.position);
			(attachment.end == 0 ? left_positions : right_positions).insert(attachment.position);
		}
		std::sort(positions.begin(), positions.end());
		positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
		vector<int> positions_with_junctions;
		for (int position : positions) {
			positions_with_junctions.push_back(position);
			if (left_positions.count(position) > 0 && right_positions.count(position) > 0) {
				positions_with_junctions.push_back(position);
			}
		}
		positions = std::move(positions_with_junctions);
		if (length == 0) {
			positions = {0, 0};
		} else {
			//links entering at the very end (or leaving at the very start) of the contig need a zero-length sub-contig there
			if (left_positions.count(length) > 0 && right_positions.count(length) == 0) {
				positions.push_back(length);
			}
			if (right_positions.count(0) > 0 && left_positions.count(0) == 0) {
				positions.insert(positions.begin(), 0);
			}
		}

		vector<string> sub_names;
		for (int sub = 0; sub + 1 < static_cast<int>(positions.size()); ++sub) {
			sub_names.push_back(sub_contig_name(contig, positions[sub], positions[sub + 1]));
			new_contig_order.push_back(sub_names.back());
			new_contig_coordinates.push_back({contig, positions[sub], positions[sub + 1]});
			new_list_of_contigs[sub_names.back()] = ContigSides{};
			if (sub > 0) {
				new_list_of_links.push_back(Link{sub_names[sub], 0, 0, sub_names[sub - 1], 1, 0});
				new_list_of_contigs[sub_names[sub]].left.push_back(static_cast<int>(new_list_of_links.size()) - 1);
			}
		}

		for (const Attachment& attachment : attachments) {
			int sub = 0;
			if (attachment.end == 0) {
				//the link enters the sub-contig starting at the breakpoint
				sub = static_cast<int>(std::lower_bound(positions.begin(), positions.end(), attachment.position) - positions.begin());
				sub = std::min(sub, static_cast<int>(sub_names.size()) - 1);
			} else {
				//the link leaves from the sub-contig ending at the breakpoint
				sub = static_cast<int>(std::upper_bound(positions.begin(), positions.end(), attachment.position) - positions.begin()) - 2;
				sub = std::max(sub, 0);
			}
			const string& new_name = sub_names[sub];

			int& new_idx = new_index_of_link[attachment.link_idx];
			if (new_idx == -1) {
				new_list_of_links.push_back(list_of_links[attachment.link_idx]);
				new_idx = static_cast<int>(new_list_of_links.size()) - 1;
				Link& new_link = new_list_of_links[new_idx];
				new_link.written = false;
				if (cut1[attachment.link_idx] != -1) {
					new_link.overlap1 = 0;
					new_link.overlap2 = 0;
				}
				if (attachment.end == 0) {
					new_list_of_contigs[new_name].left.push_back(new_idx);
				} else {
					new_list_of_contigs[new_name].right.push_back(new_idx);
				}
			}
			(attachment.side == 1 ? new_list_of_links[new_idx].name1 : new_list_of_links[new_idx].name2) = new_name;
			number_of_sides_rewired[attachment.link_idx] += 1;
		}
	}

	//checks: both sides of every link have been rewired, and the cuts of each link remove exactly its overlap (otherwise the
	//neighbour would keep an overlap that is written as 0M and sequence would be lost or duplicated)
	int number_of_links_not_blunt = 0;
	int number_of_links_badly_rewired = 0;
	int number_of_imperfect_links = 0;
	for (int link_idx = 0; link_idx < number_of_links; ++link_idx) {
		if (!link_is_valid[link_idx]) {
			continue;
		}
		const Link& link = list_of_links[link_idx];
		if (number_of_sides_rewired[link_idx] != 2) {
			number_of_links_badly_rewired += 1;
		}
		if (cut1[link_idx] == -1) {
			number_of_links_not_blunt += 1;
		} else if (!is_hairpin(link) && cut1[link_idx] + cut2[link_idx] != link.overlap1) {
			number_of_links_badly_rewired += 1;
		} else if (link.overlap1 != link.overlap2) {
			number_of_imperfect_links += 1;
		}
	}
	if (number_of_links_badly_rewired > 0) {
		std::cout << "WARNING: in bluntify, " << number_of_links_badly_rewired
				  << " links were not rewired correctly, some sequence may be lost or duplicated. This is a bug, please report it" << std::endl;
	}
	if (number_of_imperfect_links > 0) {
		std::cout << "WARNING: in bluntify, " << number_of_imperfect_links
				  << " links have an imperfect overlap, some bases may be duplicated or lost around them" << std::endl;
	}
	if (number_of_links_not_blunt > 0) {
		std::cout << "WARNING: in bluntify, " << number_of_links_not_blunt
				  << " links could not be made blunt (overlaps of short contigs overlap each other), they are kept with their overlap" << std::endl;
	}

	ofstream fo(gfa_out);
	if (!fo) {
		throw std::runtime_error("Cannot open output GFA: " + gfa_out);
	}
	ifstream open2(gfa_in);
	if (!open2) {
		throw std::runtime_error("Cannot open input GFA: " + gfa_in);
	}

	//sub-contigs of a same contig are consecutive in new_contig_order: read each original sequence only once
	string current_original_name;
	string seq;
	bool has_current_original = false;
	for (int i = 0; i < static_cast<int>(new_contig_order.size()); ++i) {
		const string& contig = new_contig_order[i];
		const SubContig& coordinates = new_contig_coordinates[i];

		if (!has_current_original || coordinates.original_name != current_original_name) {
			current_original_name = coordinates.original_name;
			has_current_original = true;
			open2.clear();
			open2.seekg(location_of_contigs_in_gfa[current_original_name]);
			string sline;
			getline(open2, sline);
			if (!sline.empty() && sline.back() == '\r') {
				sline.pop_back();
			}
			vector<string> ls = split_tab(sline);
			seq = (ls.size() > 2) ? std::move(ls[2]) : string();
		}

		if (coordinates.end - coordinates.start >= short_contig_length) {
			fo << "S\t" << contig << "\t";
			fo.write(seq.data() + coordinates.start, coordinates.end - coordinates.start);
			fo << "\t";
			auto cov = coverage.find(current_original_name);
			if (cov != coverage.end()) {
				fo << cov->second << "\n";
			} else {
				fo << "\n";
			}
		}
	}

	for (const string& contig : new_contig_order) {
		for (int link_idx : new_list_of_contigs[contig].left) {
			if (!new_list_of_links[link_idx].written) {
				write_link(fo, new_list_of_links[link_idx]);
				new_list_of_links[link_idx].written = true;
			}
		}

		for (int link_idx : new_list_of_contigs[contig].right) {
			if (!new_list_of_links[link_idx].written) {
				write_link(fo, new_list_of_links[link_idx]);
				new_list_of_links[link_idx].written = true;
			}
		}
	}
}

void remove_contigs_of_length_0(const std::string& gfa_in, const std::string& gfa_out) {
	unordered_map<string, bool> contig_is_empty;
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
				const string& name = fields[1];
				if (contig_is_empty.find(name) == contig_is_empty.end()) {
					contig_order.push_back(name);
				}
				contig_is_empty[name] = fields.size() <= 2 || fields[2].empty();
				contig_lines[name] = line + "\n";
			} else if (line.rfind("L", 0) == 0) {
				links.push_back(split_tab(line));
			}
		}
	}

	unordered_set<string> zero_length_contigs;
	for (const auto& kv : contig_is_empty) {
		if (kv.second) {
			zero_length_contigs.insert(kv.first);
		}
	}

	unordered_map<string, unordered_set<Neighbor, NeighborHash>> left_neighbors;
	unordered_map<string, unordered_set<Neighbor, NeighborHash>> right_neighbors;
	for (const auto& kv : contig_is_empty) {
		left_neighbors[kv.first] = {};
		right_neighbors[kv.first] = {};
	}

	//links that are not blunt (could not be made blunt by fancier_overlap_removal) keep their CIGAR, the others are output as 0M
	unordered_map<WrittenLink, string, WrittenLinkHash> non_blunt_cigars;

	for (const vector<string>& fields : links) {
		if (fields.size() < 5) {
			continue;
		}
		const string& from_name = fields[1];
		const string& to_name = fields[3];
		const string& from_orient = fields[2];
		const string& to_orient = fields[4];
		if (fields.size() > 5 && fields[5] != "0M" && fields[5] != "*") {
			string to_end = (to_orient == "+") ? "-" : "+";
			non_blunt_cigars[WrittenLink{from_name, from_orient, to_name, to_end}] = fields[5];
			non_blunt_cigars[WrittenLink{to_name, to_end, from_name, from_orient}] = fields[5];
		}
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

	//neighbours of each end of the non-empty contigs, going through any number of zero-length contigs (also through
	//self-links of zero-length contigs, e.g. a hairpin left at a junction by fancier_overlap_removal)
	if (!zero_length_contigs.empty()) {
		auto neighbors_through_zero_length_contigs = [&](const unordered_set<Neighbor, NeighborHash>& direct_neighbors) {
			unordered_set<Neighbor, NeighborHash> result;
			unordered_set<Neighbor, NeighborHash> visited_zero_contigs; //zero-length contig and end through which it was entered
			vector<Neighbor> stack(direct_neighbors.begin(), direct_neighbors.end());
			while (!stack.empty()) {
				Neighbor neighbor = stack.back();
				stack.pop_back();
				if (zero_length_contigs.find(neighbor.name) == zero_length_contigs.end()) {
					result.insert(neighbor);
				} else if (visited_zero_contigs.insert(neighbor).second) {
					//entered by its left end ("-"): leave by its right end, and conversely
					const unordered_set<Neighbor, NeighborHash>& next = (neighbor.orient == "-") ? right_neighbors[neighbor.name]
																								  : left_neighbors[neighbor.name];
					stack.insert(stack.end(), next.begin(), next.end());
				}
			}
			return result;
		};

		unordered_map<string, unordered_set<Neighbor, NeighborHash>> new_left_neighbors;
		unordered_map<string, unordered_set<Neighbor, NeighborHash>> new_right_neighbors;
		for (const string& name : contig_order) {
			if (zero_length_contigs.find(name) == zero_length_contigs.end()) {
				new_left_neighbors[name] = neighbors_through_zero_length_contigs(left_neighbors[name]);
				new_right_neighbors[name] = neighbors_through_zero_length_contigs(right_neighbors[name]);
			}
		}
		left_neighbors = std::move(new_left_neighbors);
		right_neighbors = std::move(new_right_neighbors);
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

		unordered_set<WrittenLink, WrittenLinkHash> written_links;
		auto cigar_of = [&](const WrittenLink& link) -> string {
			auto it = non_blunt_cigars.find(link);
			return (it == non_blunt_cigars.end()) ? string("0M") : it->second;
		};

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
					   << ((to.orient == "+") ? "-" : "+") << "\t" << cigar_of(link_tuple) << "\n";
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
					   << ((to.orient == "+") ? "-" : "+") << "\t" << cigar_of(link_tuple) << "\n";
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
 * Each end loses min(trim_isolated_contigs_length, length/2) bases; isolated contigs that are left with no sequence are dropped
 * (they have no link, so nothing else changes).
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
        int number_of_dropped_isolated_contigs = 0;
        while (getline(fi, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (line.rfind("S", 0) == 0) {
                vector<string> fields = split_tab(line);
                if (fields.size() < 2) {
                    continue;
                }
                const string& name = fields[1];
                if (isolated_contigs.find(name) != isolated_contigs.end()) {
                    if (fields.size() < 3) {
                        fields.push_back("");
                    }
                    string& seq = fields[2];
                    int trim_length = std::min(trim_isolated_contigs_length, static_cast<int>(seq.size()) / 2);
                    seq = seq.substr(static_cast<size_t>(trim_length), seq.size() - 2 * static_cast<size_t>(trim_length));
                    if (seq.empty()) {
                        //fully trimmed (even length <= 2*trim_isolated_contigs_length): never output an empty contig
                        number_of_dropped_isolated_contigs += 1;
                        continue;
                    }
                    fo << join_tab(fields) << "\n";
                } else {
                    fo << line << "\n";
                }
            } else {
                fo << line << "\n";
            }
        }
        if (number_of_dropped_isolated_contigs > 0) {
            std::cout << "In bluntify, " << number_of_dropped_isolated_contigs
                      << " short isolated contigs were entirely trimmed and dropped" << std::endl;
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

	int trim_isolated_contigs_length = 0;
	string tmpdir = ".";
	vector<string> positional;

	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];
		if (arg == "-t" || arg == "--trim_isolated") {
			if (i + 1 >= argc) {
				std::cerr << "Error: " << arg << " requires a value\n";
				return 1;
			}
			try {
				size_t parsed = 0;
				trim_isolated_contigs_length = std::stoi(argv[++i], &parsed);
				if (parsed != string(argv[i]).size() || trim_isolated_contigs_length < 0) {
					throw std::invalid_argument("negative or not an integer");
				}
			} catch (const std::exception&) {
				std::cerr << "Error: " << arg << " requires a non-negative integer, got " << argv[i] << "\n";
				return 1;
			}
		} else if (arg == "-n" || arg == "--no_overlaps") {
			//used to be passed (as a bool!) as the trimming length of isolated contigs; now ignored
			std::cerr << "Warning: " << arg << " is deprecated and ignored, use --trim_isolated LENGTH to trim isolated contigs\n";
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
				  << " input output [-t|--trim_isolated LENGTH] [--tmpdir TMPDIR] [--version]\n";
		return 1;
	}

	bluntify(positional[0], positional[1], trim_isolated_contigs_length, tmpdir);
	return 0;
}
