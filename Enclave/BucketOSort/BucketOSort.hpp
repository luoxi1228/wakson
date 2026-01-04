#pragma once
#ifndef __BUCKET_OSORT__
#define __BUCKET_OSORT__

  #ifndef BEFTS_MODE  
    #include "../../Globals.hpp"
    #include "../../CONFIG.h"
    #include "../Globals_Enclave.hpp"
    #include "../oasm_lib.h"
    #include "../utils.hpp"
    #include "ObliviousPrimitives.hpp"
    #include "MergeSplitNode.hpp"
    #include "FilterNode.hpp"
  #endif

  double DecryptAndBOS(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, 
          bos_params* params, unsigned char *result_buffer);
  
  double BOS(unsigned char *buffer, size_t N, size_t block_size, bos_params *params,
            unsigned char* &output, size_t nthread, unsigned char *result_buf, enc_ret *ret);

  double BORP(unsigned char *input, size_t N, size_t block_size, bos_params* params,
            unsigned char* &expanded_buffer, size_t nthreads, unsigned char *result_buf, enc_ret *ret);

  template<OSwap_Style oswap_style>
  void RecursiveShuffle_M2_inner_parallel(unsigned char *buf, uint64_t N, size_t block_size, bool *selected_list, threadid_t thread_id, size_t nthreads);

  #include "BucketOSort.tcc"

#endif
