#include "BOS_TC.hpp"

size_t MarkReals(unsigned char *bucket_ptr, uint64_t Z, size_t block_size, bool *selected) {
  unsigned char *ptr = bucket_ptr;
  uint64_t destination_label;
  size_t num_reals = 0;
  for(uint64_t i=0; i<Z; i++) {
    destination_label = *((uint64_t*) ptr);
    bool is_real = !(destination_label==UINT64_MAX);
    selected[i]=is_real;
    num_reals+=is_real;
    ptr+=block_size;
  }
  return num_reals; 
}

void StripLabels(unsigned char *bucket_ptr, uint64_t Z, size_t block_size_with_tag){
  unsigned char *old_ptr = bucket_ptr;
  unsigned char *new_ptr = bucket_ptr; 

  for(uint64_t i = 0; i<Z; i++) {
    memmove(new_ptr, old_ptr+8, block_size_with_tag-8); 
    new_ptr+=(block_size_with_tag-8);
    old_ptr+=block_size_with_tag;
  } 
}


void MarkForCompaction(unsigned char *buffer, int bit_position, uint64_t Z, size_t block_size, bool * selected) {
  //if (layer-th MSB is NOT set): Mark selected for it to 1
  //else Mark selected for it to 0

  unsigned char *bfr_ptr = buffer;

  uint64_t going_left = 0;
  uint64_t going_right = 0;
  uint64_t bucket_label;

  for(uint64_t i=0; i<2*Z; i++) {
    bucket_label = *((uint64_t*)bfr_ptr);
    bool direction_bit = bucket_label & (1<<(bit_position-1)); 
    // If bit at bit_posistion is set, the packet goes right.
    bool is_real = (bucket_label!= UINT64_MAX);
    going_right += (direction_bit & is_real);
    going_left += (!(direction_bit) & is_real);
    bfr_ptr+=block_size;
  }

  uint64_t dummies_to_mark = Z - going_left;
  bfr_ptr = buffer;
  for(uint64_t i=0; i<2*Z; i++) {
    bucket_label = *((uint64_t*)bfr_ptr);
    bool direction_bit = bucket_label & (1<<(bit_position-1));
    bool is_real = (bucket_label!= UINT64_MAX);

    bool mark_dummy_left = (!is_real & dummies_to_mark!=0);
    dummies_to_mark-=mark_dummy_left;
    selected[i] = (is_real & !direction_bit)|mark_dummy_left;

    bfr_ptr+=block_size;
  }
}

double BORP_TC(unsigned char *decrypted_buffer, uint64_t B, uint64_t Z, size_t block_size,
  unsigned char *result_buffer) {
  double ptime; 
  
  if(block_size == 16) {
    RouteThroughBRN<OSWAP_16X>(decrypted_buffer, B, Z, block_size);
    ProcessFinalLayer<OSWAP_8>(decrypted_buffer, B, Z, block_size, result_buffer);
  } else if (block_size%16==0) {
    RouteThroughBRN<OSWAP_16X>(decrypted_buffer, B, Z, block_size);
    ProcessFinalLayer<OSWAP_8_16X>(decrypted_buffer, B, Z, block_size, result_buffer);
  } else {
    RouteThroughBRN<OSWAP_8_16X>(decrypted_buffer, B, Z, block_size);
    ProcessFinalLayer<OSWAP_16X>(decrypted_buffer, B, Z, block_size, result_buffer);
  }

  return ptime;
}

double DecryptAndBORP_TC(unsigned char *encrypted_buffer, uint64_t N, 
        size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret) {

  long t1, t2;
  double compare_ms;
  // Using static bucket size Z = 512, for failure probability = e^-{Z/6} 
  //uint64_t Z = 512;
  uint64_t Z = 512;
  // N_prime = next 2 power >= N (BORP_TC only works for powers of 2)
  uint64_t N_prime = pow2_gt(N);
  uint64_t B = ceil(float(2 * N_prime)/Z);
  if(B==0)
    B=1;

  // Decrypt buffer to decrypted_buffer
  unsigned char *decrypted_buffer = NULL;
  size_t data_size = encrypted_block_size - 16 - 12;
  unsigned char *random_bytes = new unsigned char[8*N];
  getBulkRandomBytes(random_bytes, 8*N);

  size_t decrypted_block_size = decryptBuffer_attachRTags_addDummies(encrypted_buffer, 
          N, N_prime, B, Z, encrypted_block_size, random_bytes, &decrypted_buffer);

  if(B<=1) {
    ocall_clock(&t1);
    PRB_pool_init(1);
    RecursiveShuffle_M2(decrypted_buffer, N, decrypted_block_size);
    StripLabels(decrypted_buffer, N, decrypted_block_size);
    ocall_clock(&t2);
    compare_ms = ((double)(t2-t1))/1000.0;
    encryptBuffer(decrypted_buffer, N, decrypted_block_size-8, result_buffer);
    delete[] random_bytes;
    PRB_pool_shutdown();
    ret->OSWAP_count = OSWAP_COUNTER;
    ret->ptime = compare_ms;
    return compare_ms;
  }

  ocall_clock(&t1);
  // NOTE: We will never have decrypted_block_size==8, since attaching rTag will add 8 bytes.
  // So minimum block_size here is 16
  PRB_pool_init(1);
  BORP_TC(decrypted_buffer, B, Z, decrypted_block_size, decrypted_buffer);
  ocall_clock(&t2);

  // Encrypt buffer to result_buffer
  encryptBuffer(decrypted_buffer, N, decrypted_block_size-8, result_buffer);
  PRB_pool_shutdown();

  // CLOCKS_PER_SEC == 1000000, so CLOCKS_PER_MS == 1000
  compare_ms = ((double)(t2-t1))/1000.0;

  free(decrypted_buffer);
  delete[] random_bytes;

  ret->OSWAP_count = OSWAP_COUNTER;
  ret->ptime = compare_ms;
  return(compare_ms);
}
