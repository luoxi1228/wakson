#ifndef __WN_HPP__
#define __WN_HPP__
  
  void DecryptAndOblivWaksmanShuffle(unsigned char *encrypted_buffer, uint32_t N,
        size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret);

  void DecryptAndOblivWaksmanSort(unsigned char *encrypted_buffer, uint32_t N,
        size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret);
  
  void DecryptAndOWSS(unsigned char *encrypted_buffer, uint32_t N, size_t block_size,
    unsigned char *result_buffer, enc_ret *ret);

  void DecryptAndDJBWaksmanShuffle(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
          unsigned char *result_buffer, enc_ret *ret);

  void DecryptAndDJBWaksmanSort(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
          unsigned char *result_buffer, enc_ret *ret);
#endif
