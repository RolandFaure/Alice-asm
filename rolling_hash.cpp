#include "rolling_hash.h"
#include "robin_hood.h"

#include <iostream>

uint64_t hash_A = 0x3c8bfbb395c60474;
uint64_t hash_C = 0x3193c18562a02b4c;
uint64_t hash_G = 0x20323ed082572324;
uint64_t hash_T = 0x295549f54be24456;

// robin_hood::unordered_map<char, uint64_t> char_to_hash = {{'A', hash_A}, {'C', hash_C}, {'G', hash_G}, {'T', hash_T}};
// robin_hood::unordered_map<char, uint64_t> char_to_hash_reverse = {{'A', hash_T}, {'C', hash_G}, {'G', hash_C}, {'T', hash_A}};

inline void rol(uint64_t &x, int k){
    x = (x << k) | (x >> (64-k));
}

inline void ror(uint64_t &x, int k){
    x = (x >> k) | (x << (64-k));
}

//f(s[i+1,i+k]) = rol(f(s[i,i+k-1]),1) ^ rol(h(s[i]),k)  ^ h(s[i+k])
inline void roll_forward(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k){
    rol(hash, 1);
    rol(old_char, k);
    hash = hash ^ old_char ^ new_char;
}
//r(s[i+1,i+k]) = ror(r(s[i,i+k-1]),1) ^ ror(h(~s[i]),1) ^ rol(h(~s[i+k]),k-1)
inline void roll_reverse(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k){
    ror(hash, 1);
    ror(old_char, 1);
    rol(new_char, k-1);
    hash = hash ^ old_char ^ new_char;
}

/**
 * @brief Advance the rolling hash by one positions
 * 
 * @param foward_hash 
 * @param reverse_hash 
 * @param k 
 * @param seq 
 * @param pos_end end of the kmer on the seq
 * @param pos_begin begin of the kmer on the seq
 * @param pos_middle position of the middle of the kmer if -6666, it is not computed
 * @param hpc 
 * @return true 
 * @return false 
 */
bool roll(uint64_t &foward_hash, uint64_t &reverse_hash, int k, std::string &seq,  size_t &pos_end, long& pos_begin, long& pos_middle, bool hpc) {

    if (pos_end == seq.size()){
        return false;
    }

    if (pos_begin < 0){
        if (pos_end == 0){
            foward_hash = 0;
            reverse_hash = 0;
            pos_begin = -k;
        }

        switch (seq[pos_end])
        {
        case 'A':
            roll_forward(foward_hash, hash_A, 0, k);
            roll_reverse(reverse_hash, hash_T, 0, k);
            break;
        case 'C':
            roll_forward(foward_hash, hash_C, 0, k);
            roll_reverse(reverse_hash, hash_G, 0, k);
            break;
        case 'G':
            roll_forward(foward_hash, hash_G, 0, k);
            roll_reverse(reverse_hash, hash_C, 0, k);
            break;
        case 'T':
            roll_forward(foward_hash, hash_T, 0, k);
            roll_reverse(reverse_hash, hash_A, 0, k);
            break;
        default:
            roll_forward(foward_hash, hash_A, 0, k);
            roll_reverse(reverse_hash, hash_T, 0, k);
        }
    }
    else{
        //case switch to avoid using a dictionary
        switch (seq[pos_end])
        {
        case 'A':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_A, hash_A, k);
                roll_reverse(reverse_hash, hash_T, hash_T, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_A, hash_C, k);
                roll_reverse(reverse_hash, hash_T, hash_G, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_A, hash_G, k);
                roll_reverse(reverse_hash, hash_T, hash_C, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_A, hash_T, k);
                roll_reverse(reverse_hash, hash_T, hash_A, k);
                break;
            default:
                roll_forward(foward_hash, hash_A, hash_A, k);
                roll_reverse(reverse_hash, hash_T, hash_T, k);
            }
            break;
        case 'C':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_C, hash_A, k);
                roll_reverse(reverse_hash, hash_G, hash_T, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_C, hash_C, k);
                roll_reverse(reverse_hash, hash_G, hash_G, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_C, hash_G, k);
                roll_reverse(reverse_hash, hash_G, hash_C, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_C, hash_T, k);
                roll_reverse(reverse_hash, hash_G, hash_A, k);
                break;
            default:
                roll_forward(foward_hash, hash_C, hash_A, k);
                roll_reverse(reverse_hash, hash_G, hash_T, k);
            }
            break;
        case 'G':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_G, hash_A, k);
                roll_reverse(reverse_hash, hash_C, hash_T, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_G, hash_C, k);
                roll_reverse(reverse_hash, hash_C, hash_G, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_G, hash_G, k);
                roll_reverse(reverse_hash, hash_C, hash_C, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_G, hash_T, k);
                roll_reverse(reverse_hash, hash_C, hash_A, k);
                break;
            default:
                roll_forward(foward_hash, hash_G, hash_A, k);
                roll_reverse(reverse_hash, hash_C, hash_T, k);
            }
            break;
        case 'T':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_T, hash_A, k);
                roll_reverse(reverse_hash, hash_A, hash_T, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_T, hash_C, k);
                roll_reverse(reverse_hash, hash_A, hash_G, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_T, hash_G, k);
                roll_reverse(reverse_hash, hash_A, hash_C, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_T, hash_T, k);
                roll_reverse(reverse_hash, hash_A, hash_A, k);
                break;
            default:
                roll_forward(foward_hash, hash_T, hash_A, k);
                roll_reverse(reverse_hash, hash_A, hash_T, k);
            }
            break;
        default:
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_A, hash_A, k);
                roll_reverse(reverse_hash, hash_T, hash_T, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_A, hash_C, k);
                roll_reverse(reverse_hash, hash_T, hash_G, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_A, hash_G, k);
                roll_reverse(reverse_hash, hash_T, hash_C, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_A, hash_T, k);
                roll_reverse(reverse_hash, hash_T, hash_A, k);
                break;
            default:
                roll_forward(foward_hash, hash_A, hash_A, k);
                roll_reverse(reverse_hash, hash_T, hash_T, k);
            }
            break;
        }
    }

    pos_end++;
    while(hpc && pos_end < seq.size() && seq[pos_end] == seq[pos_end-1]){
        pos_end++;
    }

    if (pos_middle != -6666){ //special code to say that we don't want to move pos_middle

        //pos middle is in the middle of the sequence. When hpc and it falls on a stretch, take the first base of the stretch for A and C and to the last one for G and T (so that the same base is chosen in both directions)
        //first move to the end of the stretch
        while(hpc && pos_middle < seq.size() && seq[pos_middle] == seq[pos_middle+1]){
            pos_middle++;
        }
        //jumt to next stretch
        pos_middle++;
        //move to end of stretch if G or T
        while(hpc && pos_middle < seq.size() && seq[pos_middle] == seq[pos_middle+1] && (seq[pos_middle] == 'G' || seq[pos_middle] == 'T')){
            pos_middle++;
        }
    }

    pos_begin++;
    if (pos_begin < 0){
        return true;
    }
    while(hpc && pos_begin < seq.size() && pos_begin >= 1 && seq[pos_begin] == seq[pos_begin-1]){
        pos_begin++;
    }

    return true;
}

/**
 * @brief Function to perform a rolling hash of the input seq
 * 
 * @param k expected size of the kmer (equal to seq if hpc is off) Throw an error if the seq is smaller or longer than k
 * @param seq 
 * @param hpc 
 * @return uint64_t 
 */
uint64_t hash_string(int k, std::string &seq, bool hpc){
    uint64_t hash_foward = 0;
    uint64_t hash_reverse = 0;
    size_t pos_end = 0;
    long pos_begin = -k;
    int actual_k = 0;
    while (roll_f(hash_foward, k, seq, pos_end, pos_begin, hpc)){
        actual_k++;
    }
    if (actual_k != k){
        std::cerr << "ERROR in rolling_hash.cpp: Sequence  " << seq << " produces a hash of length " << actual_k << " and not " << k << std::endl;
        exit(1);
    }

    return hash_foward;
}

//same as roll but only forward hash to save time
bool roll_f(uint64_t &foward_hash, int k, std::string &seq, size_t &pos_end, long int &pos_begin, bool hpc){
        
    if (pos_end == seq.size()){
        return false;
    }

    if (pos_begin < 0){
        if (pos_end == 0){
            foward_hash = 0;
            pos_begin = -k;
        }

        switch (seq[pos_end])
        {
        case 'A':
            roll_forward(foward_hash, hash_A, 0, k);
            break;
        case 'C':
            roll_forward(foward_hash, hash_C, 0, k);
            break;
        case 'G':
            roll_forward(foward_hash, hash_G, 0, k);
            break;
        case 'T':
            roll_forward(foward_hash, hash_T, 0, k);
            break;
        default:
            std::cout << "WARNING: Non-ACGT character found in sequence, replacing them by A: "  << pos_begin<< seq << std::endl;
            roll_forward(foward_hash, hash_A, 0, k);
            break;
        }
    }
    else{
        //case switch to avoid using a dictionary
        switch (seq[pos_end])
        {
        case 'A':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_A, hash_A, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_A, hash_C, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_A, hash_G, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_A, hash_T, k);
                break;
            default:
                std::cout << "WARNING: Non-ACGT character found in sequence, replacing them by A: "  << pos_begin<< seq << std::endl;
                roll_forward(foward_hash, hash_A, hash_A, k);
                break;
            }
            break;
        case 'C':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_C, hash_A, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_C, hash_C, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_C, hash_G, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_C, hash_T, k);
                break;
            default:
                std::cout << "WARNING: Non-ACGT character found in sequence, replacing them by A: "  << pos_begin<< seq << std::endl;
                roll_forward(foward_hash, hash_C, hash_A, k);
                break;
            }
            break;
        case 'G':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_G, hash_A, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_G, hash_C, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_G, hash_G, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_G, hash_T, k);
                break;
            default:
                std::cout << "WARNING: Non-ACGT character found in sequence, replacing them by A: "  << pos_begin<< seq << std::endl;
                roll_forward(foward_hash, hash_G, hash_A, k);
                break;
                // exit(1);
            }
            break;
        case 'T':
            switch (seq[pos_begin])
            {
            case 'A':
                roll_forward(foward_hash, hash_T, hash_A, k);
                break;
            case 'C':
                roll_forward(foward_hash, hash_T, hash_C, k);
                break;
            case 'G':
                roll_forward(foward_hash, hash_T, hash_G, k);
                break;
            case 'T':
                roll_forward(foward_hash, hash_T, hash_T, k);
                break;
            default:
                std::cout << "WARNING: Non-ACGT character found in sequence, replacing them by A: "  << pos_begin<< seq << std::endl;
                roll_forward(foward_hash, hash_T, hash_A, k);
                break;
            }
            break;
        default:
            std::cout << "ERROR: Non-ACGT character found in sequence: "  << pos_begin<< seq << std::endl;
            exit(1);
        }
    }

    pos_end++;
    while(hpc && pos_end < seq.size() && seq[pos_end] == seq[pos_end-1]){
        pos_end++;
    }

    pos_begin++;
    if (pos_begin < 0){
        return true;
    }
    while(hpc && pos_begin < seq.size() && pos_begin >= 1 && seq[pos_begin] == seq[pos_begin-1]){
        pos_begin++;
    }
    
    return true;
}