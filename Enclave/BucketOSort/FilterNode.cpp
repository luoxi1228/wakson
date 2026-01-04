
#include "FilterNode.hpp"

#ifdef BOS_OUTSIDE_PRM_STORAGE
  void FilterNode::setStorageKey(unsigned char* key) {
    memcpy(storage_key, key, SGX_AESGCM_KEY_SIZE);
  }

  void FilterNode::getBucketTags(unsigned char *tags) {
    memcpy(tags, bucket_tags, f * SGX_AESGCM_MAC_SIZE);
  }
#endif

#ifdef REMOVE_FAKES_1
  int FilterNode::copyRealPackets(unsigned char *buffer){
    unsigned char *buffer_ptr = buffer;
    int total_real_packets = 0;
    for(int j = 0; j<f; j++) {
      unsigned char *bkt_bf_base_ptr = bucket_buffer + j * PARAM_Z * PACKET_SIZE;
      unsigned char *bkt_bf_ptr = bkt_bf_base_ptr; 
      bool is_packet_fake = false;
      int num_real_packets = 0;      
      int ctr = 0;

      while(is_packet_fake==false && ctr < PARAM_Z) {
        is_packet_fake=(isDummy(bkt_bf_ptr));
        if(!is_packet_fake)
          num_real_packets++;
 
        bkt_bf_ptr+=PACKET_SIZE;
        ctr+=1;
      }
      //printf("copyRealPackets: num_real_packets = %d\n", num_real_packets);

      memcpy(buffer_ptr, bkt_bf_base_ptr, num_real_packets * PACKET_SIZE); 
      buffer_ptr+=num_real_packets * PACKET_SIZE;
      total_real_packets+=num_real_packets;
    }
    //printf("real_packets_seen_so_far = %ld\n", real_packets_seen_so_far);
    return total_real_packets;
  }
#endif
