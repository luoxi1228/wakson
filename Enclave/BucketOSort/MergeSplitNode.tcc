
#ifndef __MERGESPLITNODE_TCC__
#define __MERGESPLITNODE_TCC__

  template <OSwap_Style oswap_style>
  class MergeSplitNode: public ProcessingNode {
    private:
 
      // Buffer for the regular MergeSplitNode
      //unsigned char buffer[PARAM_B * ORP_PACKET_SIZE];     
      unsigned char *buffer;
      uint64_t packets_seen_so_far;
      size_t block_size;
      size_t ORP_packet_size;
      bool multi_threaded;
      sgx_thread_mutex_t pp_mutex;

    public:
      MergeSplitNode(bos_params *params, bool multi_threaded);
      MergeSplitNode();
      ~MergeSplitNode();

      virtual void initializeBuffer(); 
      virtual int setEvictionStreams(ProcessingNode **nodes);
      virtual ProcessingNode** getEvictionStreams();
      virtual void setOutbufPtr(unsigned char *outbuf_ptr);

      virtual int getBufferOccupancy();
      /* Process a packet and then call processAndPropogate for the next node (based
          on the current_evict_stream).
      */
      virtual void processAndPropogate(unsigned char* packet_in_out);

      void process(unsigned char* packet_to_process);

      // Debug/display functions
      virtual void printNodePosition();
  };

  template<OSwap_Style oswap_style>
  MergeSplitNode<oswap_style>::MergeSplitNode() {
    evict_start_param = 0;
    packets_seen_so_far = 0;
    current_evict_stream = 1;   
    multi_threaded = false;
    initializeBuffer(); 
  }

  template<OSwap_Style oswap_style>
  MergeSplitNode<oswap_style>::MergeSplitNode(bos_params *params, bool multi_threaded) : multi_threaded(multi_threaded) {
    FOAV_SAFE_CNTXT(MSN_constructor, multi_threaded)
    f = params->f;
    B = params->B;
    d = params->d;
    evict_start_param = params->evict_start;

    log2f = params->log2f;
    node_in_layer = params->node_in_layer;
    layer = params->layer;
    ORP_packet_size = params->ORP_packet_size;
    block_size = ORP_packet_size - 16;

    packets_seen_so_far = 0;
    current_evict_stream = 0;   
   
    // Set mask for this node
    // Mask here is to extract the correct bits that this node in layer needs to
    // look at to assign it appropriate eviction stream
    mask = 0;
    //Insert log2f 1s into the LSB of mask
    for(int k=0; k<log2f; k++){
      FOAV_SAFE2_CNTXT(MSN_constructor, k, log2f)
      mask=(mask<<1)|1;
    } 
    // Now move this mask to the right log2f bits: 
    // d * log2f bits to look at in an ORP label,
    // A node at layer l needs to look at the lth log2f bits of its ORP label
    // so we shift by (d-1 - layer)*log2f bits 
    mask = mask<<(log2f*(d-1-layer));

    FOAV_SAFE_CNTXT(MSN_constructor, multi_threaded)
    if (multi_threaded) {
      sgx_thread_mutex_init(&(pp_mutex), NULL); 
    }

    buffer = (unsigned char*) malloc (B * ORP_packet_size);

    initializeBuffer(); 
  }

  template<OSwap_Style oswap_style>
  MergeSplitNode<oswap_style>::~MergeSplitNode(){
    FOAV_SAFE2_CNTXT(MSN_destructor, buffer, eviction_streams)
    if(buffer)
      free(buffer);
    if(eviction_streams)
      free(eviction_streams);
  }


  template<OSwap_Style oswap_style>
  void MergeSplitNode<oswap_style>::processAndPropogate(unsigned char* packet_to_process) {
    // Process the packet it gets in, and invokes processAndPropogate for the node
    // that is the current_evict_stream. 
    unsigned char *buffer_ptr = buffer;

    #ifdef DEBUG_PAP
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
   
    #ifdef DEBUG_PAP
      printf("Received packet: ");
      displayORPPacket(packet_to_process, block_size);
      printf("Received packet's ORP_destination_label = %ld\n", pp_destination);
      FOAV_SAFE_CNTXT(MSN_PAP, set_eviction_stream_dummy)
      if(!set_eviction_stream_dummy)
        printf("d = %d, layer = %d, log2f = %d, mask = %d, Eviction stream extracted = %ld\n", d, layer, log2f, mask, pp_eviction_stream);
    #endif

    oset_value(&pp_eviction_stream, (uint64_t)(UINT64_MAX), set_eviction_stream_dummy); 
    ((uint64_t*)(packet_to_process))[0]=pp_eviction_stream; 

    FOAV_SAFE_CNTXT(MSN_PPT, multi_threaded)
    if (multi_threaded) {
      sgx_thread_mutex_lock(&(pp_mutex));
    }

    int ret;
    //printf("packets_seen_so_far = %d, evict_start_param=%d\n", packets_seen_so_far, evict_start_param);

    FOAV_SAFE2_CNTXT(MSN_PAP, packets_seen_so_far, evict_start_param)
    if(packets_seen_so_far>=evict_start_param) {
      FOAV_SAFE_CNTXT(MSN_PAP, B)
      for(int i=0; i<B; i++) {
        FOAV_SAFE2_CNTXT(MSN_PAP, i, B)
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
        FOAV_SAFE2_CNTXT(MSN_PAP, i, B)
      }

      #ifdef DEBUG_BORP_FAILURE
        pp_eviction_stream = ((uint64_t*)(packet_to_process))[0]; 
        // NOTE: pp_eviction_stream is infact sensitive, but we use this particular piece of code
        //       only for debugging / analysing failures.
        FOAV_SAFE2_CNTXT(MSN_PAP, current_evict_stream, pp_eviction_stream)
        if(pp_eviction_stream!=-1 && current_evict_stream!=pp_eviction_stream) {
          #ifdef BEFTS_MODE 
            printf("FAIL - Incorrect eviction (most likely due to buffer overflow).\n");
          #else
            printf("REALLY? Incorrect eviction. current_evict_stream = %ld, packet's eviction stream = %ld\n", current_evict_stream, pp_eviction_stream);
          #endif
        }
      #endif

      #ifdef DISPLAY_BUFFER_OCCUPANCY 
        int buffer_occupancy = getBufferOccupancy();
        printNodePosition();
        printf("current_evict_stream = %d, buffer_occupancy = %d\n", 
                current_evict_stream, buffer_occupancy);
      #endif   
      
      ret = current_evict_stream;
      current_evict_stream = ((current_evict_stream+1)%f==0)? 0: (current_evict_stream+1);

      #ifdef DEBUG_PAP
        printf("Evicted Packet for PAP to ");
        eviction_streams[ret]->printNodePosition();
        displayORPPacket(packet_to_process, block_size);
      #endif

      // Release mutex before propogating packet to next layer    
      packets_seen_so_far+=1;
      FOAV_SAFE_CNTXT(MSN_PAP, multi_threaded)
      if (multi_threaded) {
        sgx_thread_mutex_unlock(&pp_mutex); 
      }

      eviction_streams[ret]->processAndPropogate(packet_to_process);
    } else{

      bool swap_done = false;
      for(int i=0; i<B; i++) {
        FOAV_SAFE2_CNTXT(MSN_PAP, i, B)
        // If packet_to_process is real packet move it into an empty slot in the buffer
        uint64_t bp_label = ((uint64_t*)buffer_ptr)[0]; 
        bool swap_flag = ((bp_label == UINT64_MAX) & (pp_eviction_stream != UINT64_MAX))
                          & !(swap_done);
        //printf("swap_flag = %d, i =%d\n", swap_flag, i);
        omove_buffer<oswap_style>(buffer_ptr, packet_to_process, ORP_packet_size, swap_flag); 
        swap_done=swap_done|swap_flag;
        buffer_ptr+=ORP_packet_size;
      }

      #ifdef DISPLAY_BUFFER_OCCUPANCY 
        int buffer_occupancy = getBufferOccupancy();
        printNodePosition();
        printf(" (Not evicting yet:) current_evict_stream = %d, buffer_occupancy = %d\n", 
                current_evict_stream, buffer_occupancy);
      #endif   

      packets_seen_so_far+=1;
      FOAV_SAFE_CNTXT(MSN_PAP, multi_threaded)
      if (multi_threaded) {
        sgx_thread_mutex_unlock(&pp_mutex); 
      }
    }
  }


  template <OSwap_Style oswap_style>
  void MergeSplitNode<oswap_style>::process(unsigned char* packet_to_process){
    unsigned char *buffer_ptr = buffer;
    
    #ifdef DEBUG_PAP
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
   
    #ifdef DEBUG_PAP
      printf("Received packet: ");
      displayORPPacket(packet_to_process, block_size);
      printf("Received packet's ORP_destination_label = %ld\n", pp_destination);
      if(!set_eviction_stream_dummy)
        printf("d = %d, layer = %d, log2f = %d, mask = %d, Eviction stream extracted = %ld\n", d, layer, log2f, mask, pp_eviction_stream);
    #endif

    oset_value(&pp_eviction_stream, (uint64_t)(UINT64_MAX), set_eviction_stream_dummy); 
    ((uint64_t*)(packet_to_process))[0]=pp_eviction_stream; 

    #ifdef MULTITHREADED
      int thread_ret = sgx_thread_mutex_lock(&(pp_mutex));
    #endif

    int ret;
    if(packets_seen_so_far>=evict_start_param) {
      for(int i=0; i<B; i++) {
        FOAV_SAFE2_CNTXT(MSN_process, packets_seen_so_far, evict_start_param)
        FOAV_SAFE2_CNTXT(MSN_process, i, B)
        // Swap packets, when either 1) the packet in the buffer is destined for the current_evict_stream
        // or 2) packet being processed is real, and NOT destined for the current evict stream (so you swap
        // it with a dummy packet (NOTE: The current algo will swap it with real packets as well, until it 
        // finally swaps it with a dummy packet and then will keep that.)

        // bp = buffer packet
        uint64_t bp_eviction_stream = ((uint64_t*)buffer_ptr)[0]; 
        bool swap_flag = (bp_eviction_stream == current_evict_stream)|
                         ((pp_eviction_stream != UINT64_MAX) & (pp_eviction_stream != current_evict_stream));
  
        //TODO: Compiling with this oswap_buffer templated breaks the make!!! 
        oswap_buffer<oswap_style>(packet_to_process, buffer_ptr, ORP_packet_size, swap_flag); 
        buffer_ptr+=ORP_packet_size;

        pp_eviction_stream = ((uint64_t*)(packet_to_process))[0];
      }
     
      #ifdef DISPLAY_BUFFER_OCCUPANCY 
        int buffer_occupancy = getBufferOccupancy();
        printNodePosition();
        printf("current_evict_stream = %d, buffer_occupancy = %d\n", 
                current_evict_stream, buffer_occupancy);
      #endif   
      
      ret = current_evict_stream;
      current_evict_stream = ((current_evict_stream+1)%f==0)? 0: (current_evict_stream+1);

      // Release mutex before propogating packet to next layer    
      packets_seen_so_far+=1;
      #ifdef MULTITHREADED
        sgx_thread_mutex_unlock(&pp_mutex); 
      #endif

    } else{

      bool swap_done = false;
      for(int i=0; i<B; i++) {
        FOAV_SAFE2_CNTXT(MSN_process, i, B)
        // If packet_to_process is real packet move it into an empty slot in the buffer
        uint64_t bp_label = ((uint64_t*)buffer_ptr)[0]; 
        bool swap_flag = ((bp_label == UINT64_MAX) & (pp_eviction_stream != UINT64_MAX))
                          & !(swap_done);
        //printf("swap_flag = %d, i =%d\n", swap_flag, i);
        //TODO: Port omove_buffer to <oswap_style> and uncomment the omove!
        //omove_buffer(buffer_ptr, packet_to_process, ORP_packet_size, swap_flag); 
        swap_done=swap_done|swap_flag;
        buffer_ptr+=ORP_packet_size;
      }

      #ifdef DISPLAY_BUFFER_OCCUPANCY 
        int buffer_occupancy = getBufferOccupancy();
        printNodePosition();
        printf(" (Not evicting yet:) current_evict_stream = %d, buffer_occupancy = %d\n", 
                current_evict_stream, buffer_occupancy);
      #endif   

      packets_seen_so_far+=1;
      #ifdef MULTITHREADED
        sgx_thread_mutex_unlock(&pp_mutex); 
      #endif  
    }
  
  }

  template <OSwap_Style oswap_style>
  int MergeSplitNode<oswap_style>::getBufferOccupancy(){
    unsigned char *buffer_ptr = buffer;
    int buffer_occupancy=0;
    #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
      printf("Buffer Out_addresses:" );
    #endif

    FOAV_SAFE_CNTXT(MSN_getBufferOccupancy, B)
    for(int i=0; i<B; i++) {
      FOAV_SAFE2_CNTXT(MSN_getBufferOccupancy, i, B)
      FOAV_SAFE_CNTXT(MSN_getBufferOccupancy, buffer_ptr)
      if( ((uint64_t*) buffer_ptr)[1]!=UINT64_MAX) {
        #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
          printf("(%ld,%ld)",((uint64_t*) buffer_ptr)[1],((uint64_t*) buffer_ptr)[0]);
        #endif
        buffer_occupancy++;
      }
      buffer_ptr+=ORP_packet_size;
    }

    #ifdef DISPLAY_BUFFER_OUT_ADDRESSES
      printf("\n");
    #endif

    return buffer_occupancy;
  }


  template <OSwap_Style oswap_style>
  void MergeSplitNode<oswap_style>::initializeBuffer(){
    FOAV_SAFE_CNTXT(MSN_initializeBuffer, buffer)
    FOAV_SAFE_CNTXT(MSN_initializeBuffer, ORP_packet_size)
    unsigned char dummy_packet[ORP_packet_size];
    FOAV_SAFE_CNTXT(MSN_initializeBuffer, dummy_packet)
    uint64_t dummy_val = UINT64_MAX;
    FOAV_SAFE_CNTXT(MSN_initializeBuffer, dummy_val)

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

    FOAV_SAFE_CNTXT(MSN_initializeBuffer, ORP_packet_size)
    for(int i=3*(sizeof(uint64_t)); i<ORP_packet_size; i++) {
      FOAV_SAFE_CNTXT(MSN_initializeBuffer, i)
      dummy_packet[i] = 'A';
    } 
    
    FOAV_SAFE_CNTXT(MSN_initializeBuffer, B)
    for(int i =0; i<B; i++){
      FOAV_SAFE2_CNTXT(MSN_initializeBuffer, i, B)
      memcpy(buffer_ptr, dummy_packet, ORP_packet_size);
      buffer_ptr+=ORP_packet_size; 
    }

    // Test Module to test Initialization state:
    /*  
    buffer_ptr = buffer;
    for(int i =0; i< PARAORP_packet_sizeenclave_displayPacket(buffer_ptr);
      buffer_ptr+=SERIALIZED_PACKET_SIZE;
    }
    */
  }

  template <OSwap_Style oswap_style>
  void MergeSplitNode<oswap_style>::printNodePosition() {
    printf("[%d,%d]", layer, node_in_layer); 
  }

  template <OSwap_Style oswap_style>
  void MergeSplitNode<oswap_style>::setOutbufPtr(unsigned char *outbuf_ptr) {
    //Only for FilterNode
  }


  //gets an array of f MergeSplitNode* that it uses to evict the processed packets to.
  template <OSwap_Style oswap_style>
  int MergeSplitNode<oswap_style>::setEvictionStreams(ProcessingNode **evict_nodes){
    // *nodes will be f * MergeSplitNode* pointers
    eviction_streams = evict_nodes; 
    return 1;
  }

  template <OSwap_Style oswap_style>
  ProcessingNode** MergeSplitNode<oswap_style>::getEvictionStreams(){
    return eviction_streams;
  }

#endif 
