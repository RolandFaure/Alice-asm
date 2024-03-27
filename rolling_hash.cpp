#include "rolling_hash.h"
#include "robin_hood.h"

#include <iostream>

uint64_t hash_A = 0x3c8bfbb395c60474;
uint64_t hash_C = 0x3193c18562a02b4c;
uint64_t hash_G = 0x20323ed082572324;
uint64_t hash_T = 0x295549f54be24456;

// robin_hood::unordered_map<char, uint64_t> char_to_hash = {{'A', hash_A}, {'C', hash_C}, {'G', hash_G}, {'T', hash_T}};
// robin_hood::unordered_map<char, uint64_t> char_to_hash_reverse = {{'A', hash_T}, {'C', hash_G}, {'G', hash_C}, {'T', hash_A}};

void rol(uint64_t &x, int k){
    x = (x << k) | (x >> (64-k));
}

void ror(uint64_t &x, int k){
    x = (x >> k) | (x << (64-k));
}

//f(s[i+1,i+k]) = rol(f(s[i,i+k-1]),1) ^ rol(h(s[i]),k)  ^ h(s[i+k])
void roll_forward(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k){
    rol(hash, 1);
    rol(old_char, k);
    hash = hash ^ old_char ^ new_char;
}
//r(s[i+1,i+k]) = ror(r(s[i,i+k-1]),1) ^ ror(h(~s[i]),1) ^ rol(h(~s[i+k]),k-1)
void roll_reverse(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k){
    ror(hash, 1);
    ror(old_char, 1);
    rol(new_char, k-1);
    hash = hash ^ old_char ^ new_char;
}

/**
 * @brief Advance the rolling hash by one character. Input pos is the number of chars already integrated into the hash.
 * 
 * @param foward_hash 
 * @param reverse_hash 
 * @param seq 
 * @param pos 
 * @return true 
 * @return false 
 */
bool roll(uint64_t &foward_hash, uint64_t &reverse_hash, int k, std::string &seq, size_t &pos ){
    
    if (pos == seq.size()){
        return false;
    }

    if (pos < k){
        if (pos == 0){
            foward_hash = 0;
            reverse_hash = 0;
        }

        switch (seq[pos])
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
        }
    }
    else{
        //case switch to avoid using a dictionary
        switch (seq[pos])
        {
        case 'A':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'C':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'G':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'T':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        default:
            std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
            exit(1);
        }
    }

    pos++;

    return true;
}

//same as roll but only forward hash
bool roll_f(uint64_t &foward_hash, int k, std::string &seq, size_t &pos ){
        
    if (pos == seq.size()){
        return false;
    }

    if (pos < k){
        if (pos == 0){
            foward_hash = 0;
        }

        switch (seq[pos])
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
            std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
            exit(1);
        }
    }
    else{
        //case switch to avoid using a dictionary
        switch (seq[pos])
        {
        case 'A':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'C':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'G':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        case 'T':
            switch (seq[pos-k])
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
                std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
                exit(1);
            }
            break;
        default:
            std::cout << "ERROR: Non-ACGT character found in sequence: " << seq << std::endl;
            exit(1);
        }
    }
    pos++;
    return true;
}