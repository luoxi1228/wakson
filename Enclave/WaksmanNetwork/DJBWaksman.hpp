#ifndef __DJB_WAKSMANNETWORK_HPP__
#define __DJB_WAKSMANNETWORK_HPP__

  #include "../Enclave_globals.h"
  #include <sgx_tcrypto.h>
  #include "../oasm_lib.h"
  #include "../utils.hpp"
  #include "../RecursiveShuffle/RecursiveShuffle.hpp"
  #include "../ObliviousPrimitives.hpp"
  #include "../SortingNetwork/SortingNetwork.hpp"

  class DJB_GCB_state{
    public:
      unsigned char *p;
      unsigned char *q;
      uint64_t *temp_array_p;
      uint64_t *temp_array_q;
      uint64_t *piinv;
      unsigned char *c;
      // To avoid allocating and using another temporary memory buffer, we'll reuse cp.
      // So whenever we use *cp in DJBWaksman.cpp, we recast it to the desired type.
      unsigned char *cp;

      // Per layer persistent variables
      // Replacements: f <- p, l <- q, L <- c
      //bool *f;
      size_t *F;
      //bool *l;
      //size_t *L;

      DJB_GCB_state(uint64_t N) {
        p = new unsigned char [N * 8];
        q = new unsigned char [N * 8];
        temp_array_p = new uint64_t[N];
        temp_array_q = new uint64_t[N];
        piinv = new uint64_t[N];

        int lgN = calculatelog2(N);
        // max_temp_size, max_ts = N * lg(N/2) (in interleaveControlBits)
        // cp, when being used as cp (and not a temp mem buff), needs to store N uint64_t
        // in the first level of recursion
        size_t size_multiplier;
        (lgN<=8) ? size_multiplier = 8: size_multiplier = lgN;
        cp = new unsigned char[N * size_multiplier];
        c = new unsigned char [N * size_multiplier];

        // Per layer
        //size_t nlgn = N * lgN;
        //size_t nby2lgn = N/2 * lgN;
        //f = new bool[N];
        F = new size_t[N*lgN];
        //l = new bool[N];
        //L = new size_t[2*N];
      }

      ~DJB_GCB_state() {
        delete [] p;
        delete [] q;
        delete [] temp_array_p;
        delete [] temp_array_q;
        delete [] piinv;
        delete [] c;
        delete [] cp;
    
        //delete [] f;
        delete [] F;
        //delete [] l;
        //delete [] L;
      } 
  };

  void generateControlBits(uint64_t *random_permutation, size_t N, bool *controlbits);

  template <OSwap_Style oswap_style>
  void applyPermutation(bool *control_bits, size_t C, unsigned char *buffer, size_t N,
          size_t block_size);

  template <OSwap_Style oswap_style>
  void applyPermutation(bool *control_bits_f, bool *control_bits_b, unsigned char *buffer, size_t N,
          size_t gN, size_t block_size, size_t offset);

  void DJBWaksmanShuffle(unsigned char *buffer, size_t N, size_t block_size, enc_ret *ret);

  void DecryptandDJBWaksmanShuffle(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret);

  /*
  void DJBWaksmanSort(unsigned char *buffer, size_t N, size_t block_size, unsigned char *result_buffer, enc_ret *ret); 

  void DJBWaksmanDecryptAndSort(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret); 
 */

  #include "DJBWaksman.tcc"

#endif
