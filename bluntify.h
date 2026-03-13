#ifndef BLUNTIFY_H
#define BLUNTIFY_H

#include <string>

std::string reverse_complement(const std::string& seq);

void basic_overlap_removal(const std::string& gfa_in, const std::string& gfa_out);
void fancier_overlap_removal(const std::string& gfa_in, const std::string& gfa_out,
							int short_contig_length = 0);
void remove_contigs_of_length_0(const std::string& gfa_in, const std::string& gfa_out);

void bluntify(const std::string& input, const std::string& output,
			  int trim_isolated_contigs_length = 0, const std::string& tmpdir = ".");

int bluntify_main(int argc, char** argv);

#endif
