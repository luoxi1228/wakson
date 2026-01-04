#ifndef __OSORT_RS_NS_HPP
#define __OSORT_RS_NS_HPP
 
  #include <stdio.h>
  #include <stdlib.h>
  #include "../utils.hpp"
  #include "../RecursiveShuffle/RecursiveShuffle.hpp"
 
  double DecryptAndOSortwithRSandNS(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer);

#endif
