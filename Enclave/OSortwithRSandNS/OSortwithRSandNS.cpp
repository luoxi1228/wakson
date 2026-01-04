
#include "OSortwithRSandNS.hpp"

double DecryptAndOSortwithRSandNS(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret) {

  // Decrypt buffer to decrypted_buffer
  unsigned char *decrypted_buffer = NULL;
  size_t decrypted_block_size = decryptBuffer(encrypted_buffer, N, encrypted_block_size, &decrypted_buffer);
  
  long start_time, qstart_time, stop_time;
  double ptime, qtime;

  //double rstime, qstime;
  //long stop_rs;

  ocall_clock(&start_time);

  // First Obliviously Shuffle the array
  PRB_pool_init(1);
  RecursiveShuffle_M2(decrypted_buffer, N, decrypted_block_size);
  //ocall_clock(&stop_rs);

  // Naive sort (qsort) the shuffled input
  ocall_clock(&qstart_time);
  qsort(decrypted_buffer, N, decrypted_block_size, compare);
  ocall_clock(&stop_time);

  //rstime = ((double)(stop_rs - start_time))/1000.0;
  //qstime = ((double)(stop_time - stop_rs))/1000.0;
  //printf("RS time = %f, QS time = %f\n", rstime, qstime);  

  ptime = ((double)(stop_time - start_time))/1000.0;
  qtime = ((double)(stop_time - qstart_time))/1000.0;

  // Encrypt buffer to result_buffer
  encryptBuffer(decrypted_buffer, N, decrypted_block_size, result_buffer);
  PRB_pool_shutdown();

  free(decrypted_buffer);

  ret->OSWAP_count = OSWAP_COUNTER; 
  ret->ptime = ptime;
  ret->qsort_time = qtime; 
  return (ptime);
}
