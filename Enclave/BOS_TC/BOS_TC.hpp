#ifndef __BOS_TC_HPP__
#define __BOS_TC_HPP__

#include "../Enclave_globals.h"
#include "../oasm_lib.h"
#include "../utils.hpp"
#include "TightCompaction/TightCompaction_v2.hpp"
#include "RecursiveShuffle/RecursiveShuffle.hpp"

size_t MarkReals(unsigned char *bucket_ptr, uint64_t Z, size_t block_size, bool *selected);

void StripLabels(unsigned char *bucket_ptr, uint64_t Z, size_t block_size);

void MarkForCompaction(unsigned char *buffer, int bit_position, uint64_t Z, size_t block_size, bool * selected);

double DecryptAndBORP(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer);

#include "BOS_TC.tcc"

#endif
