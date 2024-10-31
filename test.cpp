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

    //ensure that all the contigs in the output are subcontigs of the reference
    int n = 0;
    for (auto c: out_contigs){
        bool found = false;
        for (auto r: ref_contigs){
            if (r.find(c) != string::npos){
                found = true;
                break;
            }
        }
        if (!found){
            n++;
        }
    }

    if (n > 0){
        cout << "ERROR: " << n << " contigs in the output are not subcontigs of the reference" << endl;
        return false;
    }
    else{
        cout << "test passed" << endl;
    }

    return true;
}
