#ifndef ROLLING_HASH_H
#define ROLLING_HASH_H

#include <cstdint>
#include "robin_hood.h"

//rol(x,k) := x << k | x >> (64-k)
inline void rol(uint64_t &x, int k);

inline void ror(uint64_t &x, int k);

//f(s[i+1,i+k]) = rol(f(s[i,i+k-1]),1) ^ rol(h(s[i]),k)  ^ h(s[i+k])
inline void roll_forward(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k);

//r(s[i+1,i+k]) = ror(r(s[i,i+k-1]),1) ^ ror(h(~s[i]),1) ^ rol(h(~s[i+k]),k-1)
inline void roll_reverse(uint64_t &hash, uint64_t new_char, uint64_t old_char, int k);

uint64_t hash_string(int k, std::string &seq, bool hpc);

bool roll(uint64_t &foward_hash, uint64_t &reverse_hash, int k, std::string &seq,  size_t &pos_end, long& pos_begin, long& pos_middle, bool hpc);

bool roll_f(uint64_t &foward_hash, int k, std::string &seq, size_t &pos_end, long int &pos_begin, bool hpc); //same as roll but only forward hash


#endif