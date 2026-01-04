#ifndef __BN_HPP__
#define __BN_HPP__
  
  void DecryptAndOblivButterflyCompact(unsigned char *encrypted_buffer, uint32_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret);

  double testButterflyCompaction(unsigned char *buffer, size_t N, size_t block_size, bool *selected_list, enc_ret *ret);

#endif
