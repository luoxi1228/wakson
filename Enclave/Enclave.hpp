
#ifndef __OSORT_ENCLAVE__
  #define __OSORT_ENCLAVE__

  #include <string.h>
  #include <vector>
  #include <unordered_map>
  
  #include <stdarg.h>
  #include <stdlib.h>
  #include <stdio.h>
  #include <stdint.h>
  #include <math.h>
  #include <sgx_tcrypto.h>
  #include "sgx_thread.h"
  #include "sgx_trts.h"
  #include <assert.h>

  #include <openssl/ec.h>
  #include <openssl/bn.h>
  #include <openssl/rsa.h>
  #include <openssl/evp.h>
  #include <openssl/err.h>
  #include <openssl/rand.h>

  #include "../CONFIG.h"
  #include "Enclave_globals.h"
  #include "oasm_lib.h"
  #include "ObliviousPrimitives.hpp"
  #include "RecursiveShuffle/RecursiveShuffle.hpp"
  #include "HashTable/AES_CTR_HashFunction.hpp"
  #include "OSortwithRSandNS/OSortwithRSandNS.hpp"
  #include "BucketOSort/BucketOSort.hpp"
 
  static sgx_thread_mutex_t test_mutex;
  static int test_ctr=0;

  // Enclave Key Loading function
  void Enclave_loadTestKeys(unsigned char inkey[16], unsigned char outkey[16]);

  // Bucket Oblivious Sort Functions
  void BucketOSort_initialize(bos_params *params); 
  void BOS_loadTestKeys(unsigned char inkey[16], unsigned char outkey[16]);
  int BucketOSort_sort(unsigned char *data_in, unsigned char *data_out, size_t buffer_size);

  // DecryptAndXXXXX functions, for outside application to directly use our internal shuffle 
  // and compact functions with encryption, since the internal functions assume decrypted buffers.
  void DecryptAndShuffleM1(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size);
  void DecryptAndShuffleM2(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, size_t nthreads);
  void DecryptAndOPTightCompact(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size);
  void DecryptAndTightCompact(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size);
 
#endif
