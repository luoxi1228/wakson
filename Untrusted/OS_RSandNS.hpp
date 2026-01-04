
#ifndef __OS_RSandNS_HPP__
#define __OS_RSandNS_HPP__

  //Decrypt and Oblivious Sort with Recursive Shuffle followed by Naive Sort
  double DecryptAndOSortwithRSandNS(unsigned char *buf, size_t N, size_t encrypted_block_size, unsigned char *result_buf);

  double DecryptAndOSortwithRSandNS(unsigned char *buf, size_t N, size_t encrypted_block_size, unsigned char *result_buf, enc_ret *ret);


#endif
