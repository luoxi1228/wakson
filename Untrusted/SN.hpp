#ifndef __SN_HPP__
#define __SN_HPP__

void OddEvenMergeSort(unsigned char *buf, uint64_t N, size_t block_size);

double DecryptAndOddEvenMergeSort(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer);

double DecryptAndOddEvenMergeSort(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer, enc_ret *ret);

double DecryptAndBitonicSort(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer, enc_ret *ret);

double DecryptAndOddEvenMergeSortShuffle(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer, enc_ret *ret);

double DecryptAndBitonicSortShuffle(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
  unsigned char *result_buffer, enc_ret *ret);

#endif
