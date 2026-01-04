#ifndef __BOS_HPP__
#define __BOS_HPP__

  //Bucket Oblivious Sort Calls:
  double DecryptAndBOS(unsigned char *encrypted_buffer, size_t N, size_t block_size, 
          bos_params *params, size_t nthreads, unsigned char *result_buffer, enc_ret *ret);
  double DecryptAndBORP(unsigned char *encrypted_buffer, size_t N, size_t block_size, 
          bos_params *params, size_t nthreads, unsigned char *result_buffer, enc_ret *ret);

  double DecryptAndBORP_TC(unsigned char *encrypted_buffer, uint64_t N, size_t block_size,
    unsigned char *result_buffer, enc_ret *ret);

  // OLD API and MT functions:
  //int BucketOSort_sort(unsigned char* data_in, unsigned char *data_out, size_t len);
  //int BucketOSort_ORoute(unsigned char* data_in, size_t len);
  //int BucketOSort_sort_MT(unsigned char *inbuf, unsigned char *outbuf, size_t buf_len);


  // MultiThreaded BORP
  double BORP_MT(unsigned char *buf, size_t N, size_t block_size, 
          bos_params *params, size_t nthreads, unsigned char *result_buf, enc_ret *ret);

  // Optimizer Tuning functions
  int process_MSN(unsigned char *packet, bos_params *params, size_t repeat);

  // Temporary Haven for functions without appropriate headers:
  double MeasureOSWAPBuffer(unsigned char *buf, size_t N, size_t block_size);


#endif
