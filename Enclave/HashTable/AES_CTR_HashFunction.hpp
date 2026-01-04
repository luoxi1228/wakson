#ifndef __AES_CTR_HASH_FUNCTION_HEADER__
#define __AES_CTR_HASH_FUNCTION_HEADER__

#include "sgx_tcrypto.h"
#include "../utils.hpp"

/*
   AES_CTR_HashFunction for replacing the non-cryptographic hash function of std::unordered_map. 
   Currently only supports uint64_t keys and uint64_t values.
 */

class AES_CTR_HashFunction{

  sgx_aes_ctr_128bit_key_t secret_key;
  unsigned char plaintext_zeroes[SGX_AESCTR_KEY_SIZE]={0};
  size_t Hash(const uint64_t inkey) const;
  __uint128_t Hash128(const uint64_t inkey) const;

  public:
  AES_CTR_HashFunction();
  size_t operator()(const uint64_t inkey) const {
    return Hash(inkey);
  }
};

#endif
