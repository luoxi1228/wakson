#ifndef __BUTTERFLY_COMPACTION_HPP__
#define __BUTTERFLY_COMPACTION_HPP__

#ifndef BEFTS_MODE
  #include <cstddef>
  #include <cstdint>
  #include <vector>
  #include "../foav.h"
  #include "../oasm_lib.h"
  #include "../ObliviousPrimitives.hpp"
  #include "../Enclave_globals.h"
  #include "../utils.hpp"
#endif


// -------- tool ----------
static inline uint32_t nextPow2_u32(uint32_t x) {
  if (x <= 1u) return 1u;
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x + 1u;
}

static inline void set_bit(uint8_t* arr, size_t idx) {
  arr[idx >> 3] |= (uint8_t)(1u << (idx & 7));
}

static inline bool get_bit_u8(const uint8_t* arr, size_t idx) {
  return ((arr[idx >> 3] >> (idx & 7)) & 1u) != 0;
}

static inline void append_bit(std::vector<uint8_t>& ctrlBits, size_t& t, bool bit) {
  const size_t needBytes = (t + 8) >> 3;
  if (ctrlBits.size() < needBytes) {
    ctrlBits.resize(needBytes, 0);
  }
  if (bit) set_bit(ctrlBits.data(), t);
  ++t;
}

static inline void append_u8(std::vector<uint8_t>& C, uint8_t bit01) {
    C.push_back((uint8_t)(bit01 & 1u));
}


static inline uint32_t range_sum_bool(const bool* A, size_t l, size_t r) {
  uint32_t s = 0;
  for (size_t i = l; i < r; ++i) s += (A[i] ? 1u : 0u);
  return s;
}

static inline uint32_t ilog2_u32(uint32_t x) {
    return 31u - (uint32_t)__builtin_clz(x);
}

static inline bool is_pow2_u32(uint32_t x) {
    return x && ((x & (x - 1u)) == 0);
}

// largest power of two <= n  (n>=1)
static inline uint32_t pow2_le_u32(uint32_t n) {
    return 1u << ilog2_u32(n);
}

// 2power 网络的 comparator 数： N/2 * log2(N)
static inline size_t count_2pow_switches(uint32_t Npow2) {
    if (Npow2 < 2) return 0;
    return ((size_t)Npow2 >> 1) * (size_t)ilog2_u32(Npow2);
}

// 网络的 comparator 数
static size_t count_or_switches(uint32_t n) {
    if (n <= 1) return 0;
    if (n == 2) return 1;
    if (is_pow2_u32(n)) return count_2pow_switches(n);

    uint32_t n1 = pow2_le_u32(n);
    uint32_t n2 = n - n1;
    return count_or_switches(n2) + count_2pow_switches(n1) + (size_t)n2;
}


// -------- offline ----------

void generateControlBits(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits);
void generateControlBitsRec(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits);
void generateControlBitsOR(const bool* selected_list, size_t n, std::vector<uint8_t>& ctrlBits);

// -------- online ----------

template <OSwap_Style oswap_style>
void applyCompaction(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits);

template <OSwap_Style oswap_style>
void applyCompactionRec(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits);

template <OSwap_Style oswap_style>
void applyCompactionRecIter(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits);

template <OSwap_Style oswap_style>
void applyCompactionOR(unsigned char* buf, size_t n, size_t block_size, const std::vector<uint8_t>& ctrlBits);


// -------- test interface ----------
double testButterflyCompaction(unsigned char* buffer, size_t N, size_t block_size, bool* selected_list, enc_ret* ret);


#include "ButterflyCompaction.tcc"
#endif
