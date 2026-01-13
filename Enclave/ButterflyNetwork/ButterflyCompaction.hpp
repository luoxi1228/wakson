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

// -------- offline ----------

void RecGenerateControlBits(const bool* selected_list,
                           size_t n,
                           std::vector<uint8_t>& ctrlBits);

// -------- online ----------
template <OSwap_Style oswap_style>
void RecApplyCompaction(unsigned char* buf,
                        size_t n,
                        size_t block_size,
                        const std::vector<uint8_t>& ctrlBits);


double testButterflyCompaction(unsigned char* buffer,
                               size_t N,
                               size_t block_size,
                               bool* selected_list,   
                               enc_ret* ret);


#include "ButterflyCompaction.tcc"
#endif
