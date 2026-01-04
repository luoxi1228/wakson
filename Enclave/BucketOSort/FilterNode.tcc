
#ifndef __FILTER_NODE_TCC__
#define  __FILTER_NODE_TCC__

  template <OSwap_Style oswap_style>
  class FilterNode : public ProcessingNode {
    private:
      //unsigned char buffer[PARAM_B_F * ORP_PACKET_SIZE];     
      unsigned char *buffer;
      size_t ORP_packet_size;
      size_t block_size;
      size_t num_packets_in_outbuf_bucket; 

      #ifdef DEBUG_FN_TO_OUTBUF
        uint64_t real_packets_seen_so_far=0;
      #endif 

      #ifdef REMOVE_FAKES_1
        #ifndef BOS_OUTSIDE_PRM_STORAGE 
          unsigned char bucket_buffer[(uint64_t)(PARAM_F * PARAM_Z * PACKET_SIZE)];     
        #endif
      #endif
      
      unsigned char *outbuf; 
      uint64_t num_out_packets;
      uint64_t packets_seen_so_far;
      
      #ifdef BOS_OUTSIDE_PRM_STORAGE   
        unsigned char storage_key[SGX_AESGCM_KEY_SIZE];
        unsigned char bucket_tags[PARAM_F * SGX_AESGCM_MAC_SIZE];
      #endif
    
      bool multi_threaded;
      sgx_thread_mutex_t pp_mutex;

    public:

      FilterNode(bos_params *params, bool multi_threaded);
      FilterNode(); 
      ~FilterNode();
 
      #ifdef BOS_OUTSIDE_PRM_STORAGE
        void setStorageKey(unsigned char *key);
        void getBucketTags(unsigned char *tags);
      #endif     

      #ifdef REMOVE_FAKES_1
        int copyRealPackets(unsigned char *buffer_ptr);
      #endif 

      void fillBucketsWithDummies();
    
      #ifdef DEBUG_FN_TO_OUTBUF
        uint64_t getRealPacketSeenSoFar();
      #endif

      virtual void initializeBuffer(); 
      virtual int setEvictionStreams(ProcessingNode **nodes);
      virtual ProcessingNode** getEvictionStreams();
      virtual void setOutbufPtr(unsigned char *outbuf_ptr);

      virtual int getBufferOccupancy();
      /* Process a packet and then call processAndPropogate for the next node (based
          on the current_evict_stream).
      */
      virtual void processAndPropogate(unsigned char* packet_in_out);

      // Debug/display functions
      virtual void printNodePosition();
  }; 


#ifdef DEBUG_FN_TO_OUTBUF
  template <OSwap_Style oswap_style>
  uint64_t FilterNode<oswap_style>::getRealPacketSeenSoFar(){
    return real_packets_seen_so_far;
  }
#endif

template <OSwap_Style oswap_style>
int FilterNode<oswap_style>::getBufferOccupancy(){
  unsigned char *buffer_ptr = buffer;
  int buffer_occupancy=0;
  #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
    printf("Buffer Out_addresses:" );
  #endif

  FOAV_SAFE_CNTXT(FN_getRealPacketSeenSoFar, B)
  for(int i=0; i<B; i++) {
    FOAV_SAFE2_CNTXT(FN_getRealPacketSeenSoFar, i, B)
    FOAV_SAFE_CNTXT(FN_getRealPacketSeenSoFar, buffer_ptr)
    if( ((uint64_t*) buffer_ptr)[1]!=UINT64_MAX) {
      #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
        printf("(%ld,%ld)",((uint64_t*) buffer_ptr)[1],((uint64_t*) buffer_ptr)[0]);
      #endif
      buffer_occupancy++;
    }
    buffer_ptr+=ORP_packet_size;
    FOAV_SAFE2_CNTXT(FN_getRealPacketSeenSoFar, i, B)
  }

  #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
    printf("\n");
  #endif

  return buffer_occupancy;
}


template <OSwap_Style oswap_style>
void FilterNode<oswap_style>::initializeBuffer(){
  
  FOAV_SAFE_CNTXT(FN_initializeBuffer, ORP_packet_size)
  unsigned char dummy_packet[ORP_packet_size];
  uint64_t dummy_val = UINT64_MAX;

  unsigned char *ptr = dummy_packet;
  unsigned char *buffer_ptr = buffer;

  // Set eviction_stream, ORP_label, and key to all dummy_val for
  // initializing buffer with dummy packets
  memcpy(ptr, (unsigned char*) &dummy_val, sizeof(uint64_t));
  ptr+=sizeof(uint64_t);
  memcpy(ptr, (unsigned char*) &dummy_val, sizeof(uint64_t));
  ptr+=sizeof(uint64_t);
  memcpy(ptr, (unsigned char*) &dummy_val, sizeof(uint64_t));
  ptr+=sizeof(uint64_t);

  for(int i=3*(sizeof(uint64_t)); i<ORP_packet_size - 3 *sizeof(uint64_t); i++) {
    FOAV_SAFE2_CNTXT(FN_initializeBuffer, i, ORP_packet_size)
    dummy_packet[i] = 'A';
  } 
 
  for(int i =0; i< B; i++){
    FOAV_SAFE2_CNTXT(FN_initializeBuffer, i, B)
    memcpy(buffer_ptr, dummy_packet, ORP_packet_size);
    buffer_ptr+=ORP_packet_size; 
  }
  #ifdef REMOVE_FAKES_1
    unsigned char *bkt_bfr_ptr = bucket_buffer;
    for(int i =0; i< PARAM_F * PARAM_Z; i++){
      FOAV_SAFE_CNTXT(FN_initializeBuffer, i)
      memcpy(bkt_bfr_ptr, dummy_packet + 2*sizeof(uint64_t), PACKET_SIZE);
      bkt_bfr_ptr+=PACKET_SIZE; 
    }
  #endif

  // Test Module to test Initialization state:
  /*  
  buffer_ptr = buffer;
  for(int i =0; i< PARAM_B; i++){
    enclave_displayPacket(buffer_ptr);
    buffer_ptr+=SERIALIZED_PACKET_SIZE;
  }
  */
}


template <OSwap_Style oswap_style>
void FilterNode<oswap_style>::setOutbufPtr(unsigned char *outbuf_ptr) {
  outbuf = outbuf_ptr;
  num_out_packets=0;
}

template <OSwap_Style oswap_style>
int FilterNode<oswap_style>::setEvictionStreams(ProcessingNode **evict_nodes){
  // *nodes will be f * ProcessingNode* pointers
  eviction_streams = evict_nodes; 
  return 1;
}

template <OSwap_Style oswap_style>
ProcessingNode** FilterNode<oswap_style>::getEvictionStreams(){
  return eviction_streams;
}


template <OSwap_Style oswap_style>
void FilterNode<oswap_style>::printNodePosition() {
  printf("[%d,%d]", layer, node_in_layer); 
}

template <OSwap_Style oswap_style>
FilterNode<oswap_style>::FilterNode(bos_params *params, bool multi_threaded)
    : multi_threaded(multi_threaded) {
  f = params->f;
  B = params->B;
  d = params->d;
  evict_start_param = params->evict_start;
  block_size = params->block_size;
  ORP_packet_size = params->ORP_packet_size;
  node_in_layer = params->node_in_layer;
  layer = params->layer;
  packets_seen_so_far = 0;
  current_evict_stream = 0; 
  log2f = params->log2f;
  num_packets_in_outbuf_bucket = params->num_packets_in_outbuf_bucket;

  // Set mask for this node
  // Mask here is to extract the correct bits that this node in layer needs to
  // look at to assign it appropriate eviction stream
  mask = 0;
  //Insert log2f 1s into the LSB of mask
  for(int k=0; k<log2f; k++){
    FOAV_SAFE2_CNTXT(FN_constructor, k, log2f)
    mask=(mask<<1)|1;
  } 
  /* This should do NOTHING: take it off
  // Now move this mask to the right log2f bits: 
  // d * log2f bits to look at in an ORP label,
  // A node at layer l needs to look at the lth log2f bits of its ORP label
  // so we shift by (d-1 - layer)*log2f bits 
  mask = mask<<(log2f*(d-1-layer));
  */

  FOAV_SAFE_CNTXT(FN_initializeBuffer, multi_threaded)
  if (multi_threaded) {
    sgx_thread_mutex_init(&(pp_mutex), NULL); 
  }

  buffer = (unsigned char*) malloc (B * ORP_packet_size);
  initializeBuffer(); 
}

template <OSwap_Style oswap_style>
FilterNode<oswap_style>::~FilterNode(){
  FOAV_SAFE_CNTXT(FN_initializeBuffer, buffer)
  if(buffer)
    free(buffer);
}

template <OSwap_Style oswap_style>
void FilterNode<oswap_style>::processAndPropogate(unsigned char* packet_to_process) {
  unsigned char *buffer_ptr = buffer;

  #ifdef DEBUG_PAP_FN
    printf("\n\nEntered PAP for node ");
    printNodePosition();
    printf("\n");
  #endif

  // pp = packet being processed, destination here is the ORP label
  uint64_t pp_destination = ((uint64_t*)(packet_to_process))[1]; 
  // Extract the log2f bits by using the mask to assign its correct eviction stream
  uint64_t pp_eviction_stream = (pp_destination & (uint64_t) mask);
  pp_eviction_stream = pp_eviction_stream >> ((d-1-layer)*log2f);
  bool set_eviction_stream_dummy = (pp_destination==UINT64_MAX);
 
  #ifdef DEBUG_PAP_FN
    printf("FN: Received packet's ORP_destination_label = %ld\n", pp_destination);
    if(!set_eviction_stream_dummy)
      printf("d = %d, layer = %d, log2f = %d, mask = %d, Eviction stream extracted = %ld\n", d, layer, log2f, mask, pp_eviction_stream);
  #endif

  oset_value(&pp_eviction_stream, UINT64_MAX, set_eviction_stream_dummy); 
  ((uint64_t*)(packet_to_process))[0]=pp_eviction_stream; 

  FOAV_SAFE_CNTXT(FN_PAP, multi_threaded)
  if (multi_threaded) {
    int thread_ret = sgx_thread_mutex_lock(&pp_mutex);
  }

  FOAV_SAFE2_CNTXT(FN_PAP, packets_seen_so_far, evict_start_param)
  if(packets_seen_so_far>=evict_start_param) {
    FOAV_SAFE_CNTXT(FN_PAP, B)
    for(int i=0; i<B; i++) {
      FOAV_SAFE_CNTXT(FN_PAP, i)
      // Swap packets, when either 1) the packet in the buffer is destined for the current_evict_stream
      // or 2) packet being processed is real, and NOT destined for the current evict stream (so you swap
      // it with a dummy packet (NOTE: The current algo will swap it with real packets as well, until it 
      // finally swaps it with a dummy packet and then will keep that.)

      // bp = buffer packet 
      uint64_t bp_eviction_stream = ((uint64_t*)buffer_ptr)[0]; 
      bool swap_flag = (bp_eviction_stream == current_evict_stream)|
                       ((pp_eviction_stream != UINT64_MAX) & (pp_eviction_stream != current_evict_stream));

      oswap_buffer<oswap_style>(packet_to_process, buffer_ptr, ORP_packet_size, swap_flag); 
      buffer_ptr+=ORP_packet_size;

      pp_eviction_stream = ((uint64_t*)(packet_to_process))[0];
      FOAV_SAFE2_CNTXT(FN_PAP, i, B)
    }
  
    #ifdef DEBUG_BORP_FAILURE
      pp_eviction_stream = ((uint64_t*)(packet_to_process))[0]; 
      // NOTE: pp_eviction_stream is infact sensitive, but we use this particular piece of code
      //       only for debugging / analysing failures.
      FOAV_SAFE2_CNTXT(FN_PAP, pp_eviction_stream, current_evict_stream)
      if(pp_eviction_stream!=-1 && current_evict_stream!=pp_eviction_stream) {
        #ifdef BEFTS_MODE 
          printf("FAIL - Incorrect eviction (most likely due to buffer overflow).\n");
        #else
          printf("REALLY?Incorrect eviction. current_evict_stream = %ld, packet's eviction stream = %ld\n", current_evict_stream, pp_eviction_stream);
        #endif
      }
    #endif
 
    #ifdef DISPLAY_BUFFER_OCCUPANCY 
      int buffer_occupancy = getBufferOccupancy();
      printf("current_evict_stream = %d, buffer_occupancy = %d\n", 
              current_evict_stream, buffer_occupancy);
    #endif   
    
    #ifdef DEBUG_PAP_FN
      printf("Evicted Packet for PAP_FN:");
      displayORPPacket(packet_to_process, block_size);
    #endif

    int ret = current_evict_stream;
    current_evict_stream = ((current_evict_stream+1)%f==0)? 0: (current_evict_stream+1);

    #ifdef DEBUG_FN_TO_OUTBUF
      if(pp_destination!=UINT64_MAX)
        printf("Was a real packet. real_packets_seen_so_far = %ld\n", ++real_packets_seen_so_far);
      else
        printf("Was a dummy. real_packets_seen_so_far = %ld\n", real_packets_seen_so_far);
      //printf("Packet %ld of bucket %ld:\n", num_out_packets, node_in_layer*f+ret);
      //displayORPPacket(packet_to_process, block_size);
      //printf("Corresponding encrypted version of packet stored:\n");
      //displayEncryptedPacket(encrypted_packet);
      //printf("\n\n"); 
    #endif
 
    unsigned char *outbuf_ptr = outbuf + 
            (ret*(num_packets_in_outbuf_bucket) + num_out_packets)*block_size;
    memcpy(outbuf_ptr, packet_to_process + 2*sizeof(uint64_t), block_size);

    FOAV_SAFE2_CNTXT(FN_PAP, ret, f)
    if(ret==f-1){
      num_out_packets+=1;
    }

    FOAV_SAFE2_CNTXT(FN_PAP, packets_seen_so_far, multi_threaded)
    packets_seen_so_far+=1;
    if (multi_threaded) {
      sgx_thread_mutex_unlock(&pp_mutex); 
    }

  } else {
    bool swap_done = false;
    for(int i=0; i<B; i++) {
      FOAV_SAFE2_CNTXT(FN_PAP, i, B)
      // If packet_to_process is real packet move it into an empty slot in the buffer
      uint64_t bp_label = ((uint64_t*)buffer_ptr)[0]; 
      bool swap_flag = ((bp_label == UINT64_MAX) & (pp_eviction_stream != UINT64_MAX))
                        & !(swap_done);
      //printf("swap_flag = %d, i =%d\n", swap_flag, i);
      omove_buffer<oswap_style>(buffer_ptr, packet_to_process, ORP_packet_size, swap_flag); 
      swap_done=swap_done|swap_flag;
      buffer_ptr+=ORP_packet_size;
    }

    packets_seen_so_far+=1;
    FOAV_SAFE_CNTXT(FN_PAP, multi_threaded)
    if (multi_threaded) {
      sgx_thread_mutex_unlock(&pp_mutex); 
    }
  }
  
}
#endif
