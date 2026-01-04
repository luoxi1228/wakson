

#ifndef __EXPANSION_HPP__
#define __EXPANSION_HPP__

#include <sgx_tcrypto.h>
#include "../oasm_lib.h"
#include "../utils.hpp"

/*
   Expansion Class:

   Expansion can take an input array of N blocks of block_size bytes each,
   with all the real blocks stored at the start of the array. It takes in a
   destinations buffer as input as well, which contain the final destinations
   for all these real packets in the buffer.
   Expand() oblivious moves these real packets to their respective desitnations
   obliviously. 
   The destinations buffer can either be an N sized buffer, if one wishes to hide
   the number of reals in the buffer, or an r sized buffer, if the input buffer 
   only contains r real blocks, and leaking r is acceptable.
   The destinations range from 0 to N-1, with trailing 0s for dummy destination
   values when an N sized destinations array is passed.

   Initialize Expansion class templated with the T (uintXX_t) type that 
   is sufficient for the use case, i.e. bits of T should be sufficient to 
   uniquely represent each of the N block buffer that is to be compacted.

   NOTE: num_threads doesn't work atm, and Expansion is limited to a single-
   threaded execution. To support multi-threading the Expansion class would 
   have to make OCAlls to outside and then re-enter with the required number of 
   threads.
 */

template <typename T> 
class Expansion {
  T N;
  T *RS_distance;
  T EN_l;

  // Function Pointers to set according to T
  void (*oset_function)(T *dest, T src, uint8_t flag);
  void (*oincrement_function)(T *value, uint8_t flag);
  uint8_t (*ogt_set_flag_function)(T value, T value_test);

  void compute_RS_distances(T r, T *destinations);
  void process_EN(uint8_t level, unsigned char *buffer_start, size_t block_size);
  void setFunctionPointers();

  public:
  void Expand(unsigned char *buf, size_t buf_len, T r, T *destinations, size_t block_size, 
      uint8_t (*isBlockReal)(unsigned char *block));
  Expansion(T N);
      // void (*oset_fn)(T *, T, uint8_t), 
      //          void (*oincrement_fn)(T*, uint8_t), uint8_t (*ogt_set_flag_fn)(T,T));
  ~Expansion();
};


template <typename T>
void inline Expansion<T>::setFunctionPointers() {
}


template <>
void inline Expansion<uint16_t>::setFunctionPointers() {
  oset_function = &(oset_uint16_t);
  oincrement_function = &(oincrement_uint16_t);
  ogt_set_flag_function = &(ogt_set_flag_uint16_t);
}


template <>
void inline Expansion<uint32_t>::setFunctionPointers() {
  oset_function = &(oset_uint32_t);
  oincrement_function = &(oincrement_uint32_t);
  ogt_set_flag_function = &(ogt_set_flag_uint32_t);
}


template <>
void inline Expansion<uint64_t>::setFunctionPointers() {
  oset_function = &(oset_uint64_t);
  oincrement_function = &(oincrement_uint64_t);
  ogt_set_flag_function = &(ogt_set_flag_uint64_t);
}


template <typename T>
void Expansion<T>::Expand(unsigned char *buf, size_t buf_len, T r, T *destinations, size_t block_size, 
      uint8_t (*isBlockReal)(unsigned char *block)){

  EN_l = calculatelog2(N)-1;
  compute_RS_distances(r, destinations);
  /*
  for(int i=0; i<r; i++){
    printf("%d, ", RS_distance[i]);
  }
  printf("\n");
  */
  for(int l=EN_l; l>=0; l--) {
    process_EN(l, buf, block_size); 
  } 

}

template <typename T>
Expansion<T>::Expansion(T N_val) {
  N = N_val;
  RS_distance = new T [N] {0};

  // Set function pointers
  setFunctionPointers();
}

template <typename T> 
Expansion<T>::~Expansion() {
  delete [] RS_distance;
}

template <typename T> 
void Expansion<T>::compute_RS_distances(T r, T *destinations){
  for(T i = 0; i<r; i++){
    uint8_t if_real = ogt_set_flag(destinations[i], 0);
    oset_function(&(RS_distance[i]), destinations[i]-i, if_real);
  }
  
} 


// Perform the OMOVEs for input level over the Tight Compaction Network
template <typename T>
void Expansion<T>::process_EN(uint8_t level, unsigned char *bfr_ptr, size_t block_size) { 

  // Number of oblivious swaps
  T comparator_dist = pow(2,level);
  T num_oswaps = N - comparator_dist;
  T sop_index = N-1;
  T fop_index = sop_index - comparator_dist;
  //T num_sections = num_oswaps/pow(2, (level));
  //T num_oswaps_per_section = num_oswaps/num_sections; 

  // bfr_fop = bfr_first_operand_pointer, bfr_sop = bfr_second_operand_pointer
  unsigned char *bfr_fop = bfr_ptr + (fop_index * block_size);
  unsigned char *bfr_sop = bfr_ptr + (sop_index * block_size);

  /*
  for(T i = num_sections; i>0; i--) {
    for(T j = num_oswaps_per_section; j>0; j--){ 
      //Adjust bfr_fop and bfr_sop for next swap
      T move_dist = RS_distance[fop_index] & ((1 << level));
      uint8_t dist_flag = ogt_set_flag(move_dist, 0);
      uint8_t move_flag = !(isDummy(bfr_fop)) && (dist_flag);
      oswap_buffer(bfr_fop, bfr_sop, block_size, move_flag);
    
      // Obliviously if move_flag, set RS_distance[thread][fop_index] to 
      // (RS_distance[thread][sop_index]-move_dist)
      RS_distance[fop_index]-= move_dist;
      oset_function(&(RS_distance[sop_index]), RS_distance[fop_index], move_flag);
      oset_function(&(RS_distance[fop_index]), 0, move_flag);

      printf("Swaps positions (%d,%d) = %d\n", fop_index, sop_index, move_flag);
      bfr_fop-=block_size;
      bfr_sop-=block_size;
      sop_index--;
      fop_index--;
    }
    
    //Adjust bfr_fop and bfr_sop for next section
    bfr_sop-=(num_oswaps_per_section)*block_size;
    bfr_fop-=(num_oswaps_per_section)*block_size;
    sop_index-=(num_oswaps_per_section);
    fop_index-=(num_oswaps_per_section);
  }
  */
  for(T i=0; i<num_oswaps; i++){ 
    T move_dist = RS_distance[fop_index] & (1 << (level));

    uint8_t dist_flag = ogt_set_flag(move_dist, 0);
    oswap_buffer(bfr_fop, bfr_sop, block_size, dist_flag);

    // Obliviously if move_flag, set RS_distance[fop_index] to 
    // (RS_distance[sop_index]-move_dist)
    RS_distance[fop_index]-= move_dist;
    oset_function(&(RS_distance[sop_index]), RS_distance[fop_index], dist_flag);
    oset_function(&(RS_distance[fop_index]), 0, dist_flag);

    bfr_fop-=block_size;
    bfr_sop-=block_size;
    sop_index--;
    fop_index--;
  }
}
#endif
