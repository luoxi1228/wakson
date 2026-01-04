
#include "AES_CTR_HashFunction.hpp"

AES_CTR_HashFunction::AES_CTR_HashFunction(){
  getRandomBytes(secret_key, SGX_AESCTR_KEY_SIZE); 
}

size_t AES_CTR_HashFunction::Hash(const uint64_t inkey) const{
  size_t output;
  sgx_status_t ret = SGX_SUCCESS;
  // SGX_AESGCM_IV_SIZE = 12 bytes
  // Use the first 32 bits as rehash counter, or just 0's in this case since the
  // unordered map would handle that.
  unsigned char ctr[SGX_AESGCM_IV_SIZE]={0};
  unsigned char ciphertext[SGX_AESCTR_KEY_SIZE]={0};
  unsigned char zeroes[SGX_AESCTR_KEY_SIZE]={0};
  // The last 64 bits are the input key to hash
  memmove(ctr + 4, (unsigned char*) &inkey, 8);
  // The Hash function never really increments ctr bits since it always does a single 
  // block of encryption, but sgx_aes_ctr_encrypt fails if we try to pass 0 bits_to_increment.
  uint32_t ctr_bits_to_increment=1; 

  ret = sgx_aes_ctr_encrypt((const sgx_aes_ctr_128bit_key_t*) &secret_key, zeroes, SGX_AESCTR_KEY_SIZE, ctr, ctr_bits_to_increment, ciphertext);
  if(ret!=SGX_SUCCESS)
    printf("return was not SGX_SUCCESS\n");

  // Return the last 64 bits of the AES_CTR ciphertext as the output of the hash function
  memmove((unsigned char*) &output, ciphertext+8, 8);

  return output;
}