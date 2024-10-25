#include "test.h"
#include "basic_graph_manipulation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::set;


/**
 * @brief compare the result to the reference, user-provided file
 * 
 * @param out_folder 
 * @param ref_file 
 * @return true 
 * @return false 
 */
bool test_assembly(string out_file, string ref_file){

    std::ifstream in_ref(ref_file);
    if (!in_ref.is_open())
    {
        cout << "Could not open file " << ref_file << endl;
        exit(1);
    }
    set<string> ref_contigs;
    string line;
    while (std::getline(in_ref, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            ref_contigs.insert(sequence);
            ref_contigs.insert(reverse_complement(sequence));
        }
    }
    in_ref.close();

    std::ifstream in_out(out_file);
    if (!in_out.is_open())
    {
        cout << "Could not open file " << out_file << endl;
        exit(1);
    }
    set<string> out_contigs;
    while (std::getline(in_out, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            out_contigs.insert(sequence);
            ref_contigs.insert(reverse_complement(sequence));
        }
    }
    in_out.close();

    bool same = true;

    //compare the two sets
    if (ref_contigs.size() != out_contigs.size())
    {
        cout << "The number of contigs in the output is different from the number of contigs in the reference" << endl;
        cout << "The number of contigs in the reference is " << ref_contigs.size() << endl;
        cout << "The number of contigs in the output is " << out_contigs.size() << endl;
        same = false;
    }

    for (auto i : ref_contigs)
    {
        if (out_contigs.find(i) == out_contigs.end())
        {
            cout << "A contig  is in the reference but not in the output" << endl;
            same = false;
        }
    }
    for (auto i : out_contigs)
    {
        if (ref_contigs.find(i) == ref_contigs.end())
        {
            cout << "A contig  is in the output but not in the reference" << endl;
            same = false;
        }
    }

    return same;
}
