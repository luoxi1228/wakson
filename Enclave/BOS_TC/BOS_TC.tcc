#ifndef __BOS_TC_TCC__
#define __BOS_TC_TCC__

template <OSwap_Style oswap_style>
void RouteThroughBRN(unsigned char *decrypted_buffer, uint64_t B, uint64_t Z, size_t block_size) {

  uint8_t logB = calculatelog2(B); 
  bool *selected = new bool[2*Z];
  unsigned char *temp_bucket = new unsigned char[block_size * Z]; 
 
  for(int i = 0; i < logB; i++) {
 
    unsigned char * buffer_ptr = decrypted_buffer;

    // Perform TC on the current layer with adjacent buckets being TC-ed
    for(int j = 0; j < B/2; j++) {  
      MarkForCompaction(buffer_ptr, logB-i, Z, block_size, selected);
      TightCompact<oswap_style>(buffer_ptr, 2*Z, block_size, selected);
      buffer_ptr+=(2*Z*block_size); 
    }

    if(i!=(logB-1)) {
      // Move the buckets around to set up the TCs for the next layer
      for(int j = 0; j<=i; j++){
        int group_size = 1<<(2+j);
        int swap_distance = (1<<(j+1))-1;
        int num_swaps_in_group = 1<<(j);
        size_t bucket_size = Z * block_size;

        for(int k=0; k<B/group_size; k++){
          for(int l=0; l<num_swaps_in_group; l++) {
            buffer_ptr = decrypted_buffer + (k * group_size + 1) * bucket_size;
            //buffer_ptr points at the second element of group
            swapBuckets(buffer_ptr, buffer_ptr + (swap_distance * bucket_size), temp_bucket, bucket_size);
            buffer_ptr+=(2*bucket_size);
          }
        }
      }
    }

  }
  delete []selected;
  delete []temp_bucket;
}

template <OSwap_Style oswap_style>
void ProcessFinalLayer(unsigned char *decrypted_buffer, uint64_t B, uint64_t Z, size_t block_size, 
        unsigned char *results) {
  //  For each bucket:
  //    Count number of reals in bucket AND Mark all reals into selected
  //    TC the bucket
  //    RS all the reals of the bucket
  //    Copy the reals into result_buffer
  size_t bucket_size_with_tags = Z * block_size; 
  unsigned char *bkt_ptr = decrypted_buffer;
  bool *selected = new bool[Z]{};
  unsigned char *result_ptr = results;
  // Update block_size to the block_size after stripping labels.
  size_t block_size_with_tags = block_size;
  block_size = block_size-8; 
  
  for (int i=0; i<B; i++) { 
    // Mark all the reals, and return the number of reals (for the RS)
    uint64_t reals = MarkReals(bkt_ptr, Z, block_size_with_tags, selected);
    // In place strip all the destination bucket labels of blocks in the bucket
    StripLabels(bkt_ptr, Z, block_size_with_tags);
    // TightCompact the reals to the head of the bucket
    TightCompact<oswap_style>(bkt_ptr, Z, block_size, selected);
    // Shuffle the reals in the bucket
    RecursiveShuffle_M2(bkt_ptr, reals, block_size);
    //Copy the reals into results
    memmove(result_ptr, bkt_ptr, reals*block_size); 

    result_ptr+=(reals*block_size);
    bkt_ptr += bucket_size_with_tags;
  }

  delete[] selected; 
}

#endif
