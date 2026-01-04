
#ifndef __BUCKET_OSORT_TCC__
#define __BUCKET_OSORT_TCC__

  template <OSwap_Style oswap_style>
  class BucketOSort{
    private:
      // number of elements to sort
      uint64_t n;
      // n_1 = 2^x s.t. n_1 < n < 2^x+1
      uint64_t n_1;
      // number of Buckets = 2n/Z
      // The b provided will be a mutliple of the f
      uint64_t b;
      // number of total (real + dummy) elements in a bucket
      uint64_t Z;
      // fan-out of the BRN
      uint64_t f;
      // depth of the BRN
      uint64_t d;
      // Buffer size for the MergeSplitNodess ( non-Filter ones)
      uint64_t B; 
      unsigned int log2f;   
      size_t block_size;
      size_t evict_start;
      size_t num_packets_in_outbuf_bucket;
      size_t ORP_packet_size;
      bos_params input_params;
      int num_threads;     
 
  
      unsigned char *ORP_dummy_packet;
      unsigned char *dummy_packet_data;

      // Needed in the MT version so that BOS class knows where to write the final output to
      unsigned char *data_out_loc; 

      #ifdef BOS_OUTSIDE_PRM_STORAGE
        unsigned char *outbuf;
        unsigned char bucket_tags[PARAM_b * SGX_AESGCM_MAC_SIZE];
        unsigned char storage_key[SGX_AESGCM_KEY_SIZE];
        unsigned char bucket[ENCRYPTED_BOS_OUTBUF_BUCKET_SIZE];
      #else
        unsigned char *outbuf;
        #ifdef MULTITHREADED
          uint32_t *reals_in_bucket;
          bool *selected_in_bucket; 
        #endif
      #endif 
      
      // Points to where the next real packet can be inserted into in outbuf
      // when performing removeFakes()
      unsigned char *rp_end_outbuf_ptr;

 
      ProcessingNode ***nodes;
   
      unsigned char decryption_key[SGX_AESGCM_KEY_SIZE];
      unsigned char encryption_key[SGX_AESGCM_KEY_SIZE];
 
      /* Gets 2 memory arrays inbuf and outbuf allocated by the untrusted
         to store the (expanded and) generated buckets with real and dummy packets
         and to store the final output of the randomized real and dummy packets at the end of
         the ORP phase */
      int getOutsidePRMMemory(unsigned char* inbuf, unsigned char* outbuf, uint64_t mem_size);
    
      /* Populate nodes with the number of MergeSplit nodes needed = d * (B/f)
         Each node needs to know which MSNode its outputs connect to (It can be agnostic of
         its input nodes?)
      */

      uint64_t calculateOutbufSize();
      void calculatelog2f();
      void generateORPDummy(unsigned char * ptr_to_packet);
  
      void flushBuffers();
      void removeFakes();
      void removeFakes_TC();
      #ifdef BOS_OUTSIDE_PRM_STORAGE
        void removeFakes_OPRM();
        void merge_OPRM(unsigned char *data, int mn_index, uint64_t pkts_per_bkt, unsigned char *temp_buff, bool isLastLevel, unsigned char *data_out);
        void mergeSort_OPRM(unsigned char *data, size_t data_size, size_t start_index, size_t end_index, unsigned char *data_out);
        void checkBucketTags(uint64_t rp_per_bkt);
      #endif

      void processPacketsThroughBRN(unsigned char *data_in); 
      void setSortParams(bos_params *params, int layer, int node_in_layer, unsigned int log2f);
      void encryptResult(unsigned char *buffer_to_write_encrypted_data);

      // Debug Functions
      void displayBuckets();
      void displayConfiguration();
      void displayOutbuf();
      void displayReals();
      void displayEncryptedReals(int rp_per_bkt);
      void displayPackets(unsigned char *buffer, uint32_t buffer_type, size_t buffer_len);
      void displaySortedData(unsigned char *data_out);
 
      static void* processPacketsThroughBRN_parallel_launch(void *voidargs);
      static void* removeFakes_p1_parallel_launch(void *voidargs);
      static void* removeFakes_p2_parallel_launch(void *voidargs);
    public:

      BucketOSort(bos_params *params, unsigned char *output_buffer, size_t nthreads);
      ~BucketOSort();
      int initializeSort();      
      //int sort(unsigned char* data_in, unsigned char* data_out, size_t len); 
      double sort(unsigned char *input, size_t N, size_t block_size, unsigned char *output, enc_ret *ret);
      double ORP(unsigned char *input, size_t N, size_t block_size, unsigned char *result_buf, enc_ret *ret);

      int ORoute(unsigned char* data_in, size_t len); 
      
      int initializeBRN();    

      // MT functions
      int setup(unsigned char *data_out);
      int remainderOS(); 
      void processPacketsThroughBRN_parallel(unsigned char *input, size_t packets_start, size_t num_packets, size_t nthreads);
      void removeFakes_p1_parallel(int bucket_start, int bucket_stop, size_t nthreads);
      void removeFakes_p2_parallel(int bucket_start, int bucket_stop, unsigned char *result_buf, size_t nthreads);
      void processPacketsThroughBRN_MT(unsigned char *inbuf_ptr, size_t packets_start, size_t num_packets); 
      void flushBuffers_MT(int depth, int msn_start, int msn_stop);
      void removeFakes_p1(int bucket_start, int bucket_stop);
      void removeFakes_p2(int bucket_start, int bucket_stop, unsigned char *result_buf, size_t nthreads);
      void TC_removeFakes_MT_p2(int bucket_start, int bucket_stop);
  }; 


  /*
    Assume initializeOutbuf() (from BucketOSort.cpp) was invoked before this with params, and the
    spawned buffer is passed in as output_buffer.
  */

  template <OSwap_Style oswap_style>
  BucketOSort<oswap_style>::BucketOSort(bos_params *params, unsigned char *output_buffer, size_t nthreads){
    n = params->n;
    n_1 = params->n_1;
    b = params->b;
    Z = params->Z;
    f = params->f;
    d = params->d;
    B = params->B;
    evict_start = params->evict_start;
    block_size = params->block_size;
    num_packets_in_outbuf_bucket = params->num_packets_in_outbuf_bucket;    
    //printf("num_packets_in_outbuf_bucket = %ld\n", num_packets_in_outbuf_bucket);
    num_threads = nthreads;
 
    memcpy(&input_params, params, sizeof(bos_params)); 

    calculatelog2f();
    input_params.log2f = log2f;

    ORP_dummy_packet = (unsigned char *) malloc (block_size + 2 *sizeof(uint64_t));
    dummy_packet_data = ORP_dummy_packet + 2*sizeof(uint64_t);
    // Setup the default dummy_packet_data for easy generation of ORP dummy packets
    for(uint64_t i =0; i<block_size-1; i++){
      FOAV_SAFE2_CNTXT(BOS_constructor, i, block_size)
      dummy_packet_data[i]='A';
    }
    dummy_packet_data[block_size-1]='\0';
    ((uint64_t*)(dummy_packet_data))[0]=UINT64_MAX;  

    setORPDummy(ORP_dummy_packet);
 
    outbuf = output_buffer;

    #ifdef MULTITHREADED
      reals_in_bucket = new uint32_t[b];
      selected_in_bucket = new bool[num_packets_in_outbuf_bucket * b];
    #endif
  }

  template <OSwap_Style oswap_style>
  BucketOSort<oswap_style>::~BucketOSort() { 
    FOAV_SAFE_CNTXT(BOS_destructor, ORP_dummy_packet)
    free(ORP_dummy_packet);
    int i = 0, j =0;
    // Go through all PN and free their buffers
    FOAV_SAFE2_CNTXT(BOS_destructor, d, b)
    FOAV_SAFE2_CNTXT(BOS_destructor, i, j)
    FOAV_SAFE_CNTXT(BOS_destructor, f)
    for(i=0; i<d; i++){
      FOAV_SAFE2_CNTXT(BOS_destructor, i, d)
      for(j=0; j<b/f; j++){
        FOAV_SAFE2_CNTXT(BOS_destructor, j, b)
        FOAV_SAFE_CNTXT(BOS_destructor, f)
        FOAV_SAFE2_CNTXT(BOS_destructor, i, d)
        if(i==d-1){
          FOAV_SAFE_CNTXT(BOS_destructor, nodes)
          FOAV_SAFE2_CNTXT(BOS_destructor, i, j)
          delete((FilterNode<oswap_style>*) nodes[i][j]);
        }
        else{
          FOAV_SAFE2_CNTXT(BOS_destructor, nodes, d)
          FOAV_SAFE2_CNTXT(BOS_destructor, i, j)
          FOAV_SAFE2_CNTXT(BOS_destructor, d, b)
          FOAV_SAFE_CNTXT(BOS_destructor, f)
          delete((MergeSplitNode<oswap_style>*) nodes[i][j]);
          FOAV_SAFE2_CNTXT(BOS_destructor, nodes, d)
          FOAV_SAFE2_CNTXT(BOS_destructor, i, j)
        }
        FOAV_SAFE2_CNTXT(BOS_destructor, j, b)
        FOAV_SAFE2_CNTXT(BOS_destructor, i, f)
      } 
      FOAV_SAFE2_CNTXT(BOS_destructor, i, d)
    }
    FOAV_SAFE_CNTXT(BOS_destructor, nodes)
    delete(nodes);
    #ifdef MULTITHREADED
      delete[] reals_in_bucket;
      delete[] selected_in_bucket;
    #endif
  }
  
 /*
  template <OSwap_Style oswap_style>
  BucketOSort<oswap_style>::~BucketOSort() { 
    FOAV_SAFE_CNTXT(BOS_destructor, ORP_dummy_packet)
    free(ORP_dummy_packet);
    int i = 0, j =0;
    //FOAV_SAFE2_CNTXT(BOS_destructor, i, d)
    // Go through all PN and free their buffers
    for(i=0; i<d; i++){
      for(j=0; j<b/f; j++){
        if(i==d-1){
          delete((FilterNode<oswap_style>*) nodes[i][j]);
        }
        else{
          delete((MergeSplitNode<oswap_style>*) nodes[i][j]);
        }
      } 
    }
    delete(nodes);
    #ifdef MULTITHREADED
      delete[] reals_in_bucket;
      delete[] selected_in_bucket;
    #endif
  }
  */
  template <OSwap_Style oswap_style>
  int BucketOSort<oswap_style>::initializeBRN() {

    /*  The BRN setup will route at Layer 0 with the MSB (log2f of them) 
        and then next layer with the next log2f MSB and so on, till the
        last layer (filter nodes) which will write them based on the log2f LSB
        into the f output buckets for that node.
    */

    // Spawn the d*(b/f) MergeSplitNodes. NOTE: for b buckets in each layer, MSN required = b/f.
    // Level = 0 to d-1, nodes in level = 0 to (b/f)-1 
    // Set level d-1 to be Filter Nodes when spawning

    // Creating mask base
    // log2f bits masked
    // Start mask with all 1s
    unsigned int mask=UINT_MAX;
    ORP_packet_size = block_size + 2 *sizeof(uint64_t);


    // Insert the log2f 0's of the mask
    FOAV_SAFE2_CNTXT(initializeBRN, log2f, mask)
    for(int l=0; l<log2f; l++){
      FOAV_SAFE2_CNTXT(initializeBRN, l, log2f)
      mask=mask<<1;
    } 

    nodes = new ProcessingNode**[d];
      
    //printf("d = %d, log2f = %d \n", d, log2f);

    for(int i=d-1; i>=0; i--) {
      // Produce mask for this layer by left shifting the masked 0 bits by log2f bits
      // Shift the 0's to the right slot by inserting as many LSB 1's required
      // The mask is originally configured correctly for d-2, (d-1) has no mask 
      FOAV_SAFE2_CNTXT(initializeBRN, i, d)
      if(i<d-2){
        FOAV_SAFE2_CNTXT(initializeBRN, i, d)
        FOAV_SAFE_CNTXT(initializeBRN, log2f)
        for(int k=0; k<log2f; k++){
          FOAV_SAFE2_CNTXT(initializeBRN, k, log2f)
          mask=(mask<<1)|1;
        } 
      }
  
      //NOTE: b%f has to be 0
      nodes[i] = new ProcessingNode*[(b/f)];
      for(unsigned int j=0; j<b/f; j++){
        FOAV_SAFE2_CNTXT(initializeBRN, b, f)
        FOAV_SAFE2_CNTXT(initializeBRN, j, d)
        if(i!=d-1){
          FOAV_SAFE2_CNTXT(initializeBRN, i, d)
          // num_packets_in_outbuf_bucket is meant for just FilterNodes, so we ignore it and set to 0 
          // for MergeSplitNodes
          setSortParams(&input_params, i, j, log2f);

          
          FOAV_SAFE_CNTXT(init_BRN, num_threads)
          nodes[i][j] = new MergeSplitNode<oswap_style>(&input_params, num_threads > 1);
          // Create the array of MergeSplitNode* of size f, which will be the next
          // nodes that this currently spawned node will be connected to in the BRN
          ProcessingNode **evictstreams = (ProcessingNode **) malloc(f*sizeof(MergeSplitNode<oswap_style>*));
          ProcessingNode **evictstreams_ptr = evictstreams;

          #ifdef PRINT_BRN_CONFIGURATION
            printf("Node[%d][%d] connects to : ",i,j);
          #endif
          
          int masked_j = mask&j;

          for(int k=0; k<f; k++){  
            FOAV_SAFE2_CNTXT(initializeBRN, k, f)
            //Get the f values of j_1,...,j_f
            // d-2 since we don't shift at all for i=d-2, ie. the LSB log2f bits is what we permute
            int m = (masked_j)|(k<<(((d-2)-i)*log2f));

            ProcessingNode *msn_ptr = nodes[i+1][m];

            memcpy(evictstreams_ptr, (unsigned char*) (&msn_ptr), sizeof(MergeSplitNode<oswap_style>*));
            evictstreams_ptr+=1;
            //evictstreams_ptr+=(sizeof(MergeSplitNode*));

            #ifdef PRINT_BRN_CONFIGURATION
              printf("[%d][%d], ",i+1, m);
            #endif
          } 
          #ifdef PRINT_BRN_CONFIGURATION
            printf("\n");
          #endif

          /*
          //Test for setting and retrieving eviction streams:       
          printf("Eviction streams set = "); 
          for(int p = 0; p<f; p++){
            printf("%ld, ", ((uint64_t*)(evictstreams))[p]);
          } 
          printf("\n");
          ProcessingNode **retrievedES;
          */
          nodes[i][j]->setEvictionStreams(evictstreams);
          /*     
          retrievedES = nodes[i][j]->getEvictionStreams();
          printf("Eviction streams retrieved = "); 
          for(int p = 0; p<f; p++){
            printf("%ld, ", ((uint64_t*)retrievedES)[p]);
          } 
          printf("\n");
          */

        } else{
          setSortParams(&input_params, i, j, log2f);
          nodes[i][j] = new FilterNode<oswap_style>(&input_params, num_threads > 1);
       
          // Tell the FilterNodes which portion of the outbuf belongs to them
          // by passing the appropiate ptr to outbuf to them
          unsigned char *outbuf_ptr = outbuf + (j*f*(num_packets_in_outbuf_bucket*block_size));
          nodes[i][j]->setOutbufPtr(outbuf_ptr);

        }
      }
    }

    return 0;
  }

  /*
    void generateORPDummy(unsigned char *ptr_to_packet) :
    Creates an ORP dummy packet ready to be encrypted. 
    It leaves the space for AES_GCM's IV and tag untouched. 
  */

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::generateORPDummy(unsigned char *ptr_to_packet) {
    // Set the ORP eviction_stream and label to UINT64_MAX  
    ((uint64_t*) ptr_to_packet)[0] = UINT64_MAX;
    ((uint64_t*) ptr_to_packet)[1] = UINT64_MAX;
    
    // Dummy packet data contains UINT64_MAX key and block_size - 8 'A's
    // Set the data to dummy data
    memcpy(ptr_to_packet + 2*sizeof(uint64_t), dummy_packet_data, block_size); 
  }


  /* Go through outbuf, and extract all the real packets and put them into
     outbuf itself */
  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes_TC() {
    //printf("In removeFakes_TC!\n");
    uint64_t real_packets=0;
    unsigned char *outbuf_ptr = outbuf;
    unsigned char *rp_end_outbuf_ptr = outbuf;

    FOAV_SAFE_CNTXT(BOS_removeFakes_TC, outbuf_ptr)
    FOAV_SAFE_CNTXT(BOS_removeFakes_TC, block_size)
    FOAV_SAFE_CNTXT(BOS_removeFakes_TC, num_packets_in_outbuf_bucket)
    unsigned char rp_current[block_size];

    bool selected[num_packets_in_outbuf_bucket]={};
    size_t reals_in_bucket = 0;
    
    //printf("In removeFakes_TC: before for loop of i over b\n");
    FOAV_SAFE_CNTXT(initializeBRN, b)
    for(uint64_t i=0; i<b; i++){
      FOAV_SAFE2_CNTXT(initializeBRN, i, b)
      reals_in_bucket=0;

      unsigned char *outbuf_ptr_scan = outbuf_ptr;
      //printf("In removeFakes_TC: before for loop of j over num_packets_in_outbuf_bucket\n");
      for(uint64_t j=0; j<num_packets_in_outbuf_bucket; j++){
        FOAV_SAFE2_CNTXT(removeFakes_TC, j, num_packets_in_outbuf_bucket)
        bool is_real = !(isDummy(outbuf_ptr_scan));
        selected[j]=is_real;
        reals_in_bucket+=is_real;
        outbuf_ptr_scan+=block_size;
      }

      #ifdef DEBUG_RFTC
        size_t count_reals = 0;
        for(uint64_t j=0; j<num_packets_in_outbuf_bucket; j++){
          FOAV_SAFE2_CNTXT(removeFakes_TC, j, num_packets_in_outbuf_bucket)
          if(selected[j]==1) {
            count_reals+=1;
          }
        }
        printf("In bucket %ld, reals = %ld\n", i, count_reals);
      #endif

      //ocall_clock(&RF_stop); 
      //RF_time = (double)(RF_stop - PAP_stop)/1000.0;
      //printf("Actual execution: PAP_time = %f, RF_time = %f\n", PAP_time, RF_time);

      // Invoke TightCompaction on each bucket
      TightCompact_v2(outbuf_ptr, num_packets_in_outbuf_bucket, block_size, selected);

      // RecursiveShuffle the realpackets to reorder them from initial order
      RecursiveShuffle_M2(outbuf_ptr, reals_in_bucket, block_size);

      memmove(rp_end_outbuf_ptr, outbuf_ptr, reals_in_bucket * block_size);
      
      outbuf_ptr = outbuf_ptr+ (block_size * num_packets_in_outbuf_bucket);
      rp_end_outbuf_ptr += (block_size * reals_in_bucket);    
      FOAV_SAFE2_CNTXT(initializeBRN, i, b)
    }
  }


  /* Go through outbuf, and extract all the real packets and put them into
     outbuf itself */
  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes() {
    uint64_t real_packets=0;
    unsigned char *outbuf_ptr = outbuf;
    unsigned char *rp_end_outbuf_ptr = outbuf;
    unsigned char rp_current[block_size];

    for(uint64_t i=0; i<b; i++){
      FOAV_SAFE2_CNTXT(removeFakes, i, b)
      bool bucket_still_has_reals = true;
      uint64_t rp_extracted_from_this_bucket = 0;

      //We only need to scan a outbuf bucket till there are no more real packets left
      while(bucket_still_has_reals){
        FOAV_SAFE_CNTXT(removeFakes, bucket_still_has_reals)
        uint32_t rp_seen_so_far=1;
        uint32_t random_coin;

        outbuf_ptr = outbuf + i * num_packets_in_outbuf_bucket * block_size + rp_extracted_from_this_bucket* block_size;    
        bool extracted_real_in_pass = false;
         
        // Set the packet we swap with to be a dummy packet
        setDummy(rp_current);

        for(uint64_t j=rp_extracted_from_this_bucket; j<num_packets_in_outbuf_bucket; j++) {
          FOAV_SAFE2_CNTXT(removeFakes, j, num_packets_in_outbuf_bucket)
          FOAV_SAFE2_CNTXT(removeFakes, j, rp_extracted_from_this_bucket)

          bool is_real =  !isDummy(outbuf_ptr);

          // make_swap with probability 1/rp_seen_so_far
          //TODO: Fix with new PRB_pool
          //getRandomBytes((unsigned char*) &random_coin, sizeof(uint32_t)); 
          uint32_t randomness = random_coin % rp_seen_so_far;
          bool make_swap = (is_real) & (randomness==0);
          oswap_buffer(rp_current, outbuf_ptr, block_size, make_swap);   

          //printf("bucket_still_has_reals = %d, extracted_real_in_pass = %d\n",
          //      bucket_still_has_reals, extracted_real_in_pass);
          outbuf_ptr+=block_size;  
          rp_seen_so_far+=is_real;
       }

        FOAV_SAFE_CNTXT(removeFakes, rp_seen_so_far)
        if(rp_seen_so_far!=1){
          memcpy(rp_end_outbuf_ptr, rp_current, block_size);
          rp_end_outbuf_ptr+=block_size;
          real_packets++;
          rp_extracted_from_this_bucket++;

          //This saves an unnecesary additional scan of the entire bucket
          if(rp_seen_so_far==2){
            bucket_still_has_reals=false;
          }
        } else{
          bucket_still_has_reals=false;
        }
      }
      #ifdef DEBUG_REMOVEFAKES
        printf("Real packets after bucket %ld = %ld\n", i, real_packets);
      #endif
    }
  }

  /*
    Send each node starting from layer 0, B*f dummy puckets
    to process to flush out their buffer of any and all real packets.
  */
  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::flushBuffers()  {
    for(unsigned int i=0; i<=d-1;i++){
      //printf("Flushing buffers in level %d\n", i);
      FOAV_SAFE2_CNTXT(flushBuffers, f, b)
      FOAV_SAFE_CNTXT(flushBuffers, i)
      for(unsigned int j=0; j<b/f; j++) {
        FOAV_SAFE2_CNTXT(flushBuffers, f, B)
        FOAV_SAFE2_CNTXT(flushBuffers, i, j)
        for(unsigned int k=0; k< B*f; k++)  {

          FOAV_SAFE2_CNTXT(flushBuffers, d, b)
          FOAV_SAFE2_CNTXT(flushBuffers, f, B)
          FOAV_SAFE2_CNTXT(flushBuffers, i, j)
          FOAV_SAFE_CNTXT(flushBuffers, k)
          unsigned char dummy_packet[ORP_packet_size]; 
          setORPDummy(dummy_packet);
          FOAV_SAFE2_CNTXT(flushBuffers, ORP_packet_size, dummy_packet)
          FOAV_SAFE_CNTXT(flushBuffers, nodes)
          nodes[i][j]->processAndPropogate(dummy_packet);
        }
      }
    }
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::setSortParams(bos_params *params, int layer, int node_in_layer, unsigned int log2f) {
    params->log2f = log2f;
    params->layer = layer;
    params->node_in_layer = node_in_layer;
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::processPacketsThroughBRN(unsigned char *data_in){
    uint64_t reals_per_bucket = n/b;
    uint64_t buckets_with_extra_real = n%b;

    // encchunk_ptr points at the encrypted data blocks passed in via data_in 
    const unsigned char *packet_ptr = data_in;
    uint64_t real_packets_processed = 0;

    setORPDummy(ORP_dummy_packet);

    //uint64_t *bucket_distribution = new uint64_t[b]{};
    FOAV_SAFE2_CNTXT(PPTBRN, b, f)
    for(unsigned int i=0; i<b/f; i++) {
      FOAV_SAFE2_CNTXT(PPTBRN, i, f)
      for(unsigned int j=0; j<f; j++){
        FOAV_SAFE2_CNTXT(PPTBRN, j, f)
        for(unsigned int k=0; k<reals_per_bucket+1; k++){
          FOAV_SAFE2_CNTXT(PPTBRN, k, reals_per_bucket)
          // New ORP Packet to be generated
          unsigned char new_ORP_packet[ORP_packet_size];
          // Reset dummy packet (in case it got overwritten with other packets in the previous 
          // processAndPropogate pass
          setORPDummy(ORP_dummy_packet);
          unsigned char *new_packet_ptr = new_ORP_packet;

          
          FOAV_SAFE2_CNTXT(PPTBRN, k, buckets_with_extra_real)
          FOAV_SAFE2_CNTXT(PPTBRN, reals_per_bucket, buckets_with_extra_real)
          FOAV_SAFE2_CNTXT(PPTBRN, i, f)
          FOAV_SAFE_CNTXT(PPTBRN, j)
          if(k<reals_per_bucket || ((k==reals_per_bucket) && (i*f+j<buckets_with_extra_real))){
            FOAV_SAFE2_CNTXT(PPTBRN, k, buckets_with_extra_real)
            FOAV_SAFE2_CNTXT(PPTBRN, reals_per_bucket, buckets_with_extra_real)
            FOAV_SAFE2_CNTXT(PPTBRN, i, f)
            FOAV_SAFE_CNTXT(PPTBRN, j)
            // Process a real packet
            //printf("\n\nProcessing packet %d\n\n", ++ctr);

            // Give the packet a random label for ORP
            uint16_t random_label;
            //TODO: Fix with new PRB_pool
            //getRandomBytes((unsigned char*) &random_label, sizeof(uint16_t));
            random_label = random_label % b;
            //bucket_distribution[random_label]++;
            ((uint64_t*) new_packet_ptr)[1] = random_label;
            // No need to set eviction_stream for ORP packet, that will be assigned by the MSN 
            // that processes this packet 
            // Shift new_packet_ptr to account for the 2 added fields for ORP packets
            new_packet_ptr = new_packet_ptr+(2*sizeof(uint64_t));
            memcpy(new_packet_ptr, packet_ptr, block_size); 

            //printf("ProcessPackets invoked with real packet label = %d\n", random_label);

            nodes[0][i]->processAndPropogate(new_ORP_packet);
            real_packets_processed=real_packets_processed+1;

            // Process a dummy packet
            nodes[0][i]->processAndPropogate(ORP_dummy_packet);
            packet_ptr+=block_size;

          }
          FOAV_SAFE2_CNTXT(PPTBRN, k, reals_per_bucket)
        }
        FOAV_SAFE2_CNTXT(PPTBRN, j, f)
      }
      FOAV_SAFE2_CNTXT(PPTBRN, i, f)
    }
    /*
    printf("b = %d\n", b);
    for (uint64_t i=0; i<b; i++){
      printf("bucket %ld got %ld real packets\n", i, bucket_distribution[i]);
    }
    delete []bucket_distribution;
    */
  }


  template <OSwap_Style oswap_style>
  double BucketOSort<oswap_style>::sort(unsigned char *input, size_t N, size_t block_size, unsigned char *output, enc_ret *ret) {

    // For now we setup a BOS object -> handled by BOS_initialize, which gets all the parameters from
    // the user/application!
    /*
    unsigned char *input_ptr=input;
    size_t block_label;
    for(size_t t=0; t<N; t++){
      memcpy((unsigned char *) &block_label, input_ptr, 8);
      printf("%ld, ",block_label);
      input_ptr+=block_size;
    }
    */

    long total_start, total_stop;
    double ttime;
    ocall_clock(&total_start);

    long PAP_stop, RF_stop;
    double PAP_time, RF_time;

    // For now this is handled by BOS_initialize call!
    initializeBRN(); 
     
    processPacketsThroughBRN(input);

    flushBuffers();
    ocall_clock(&PAP_stop);
    PAP_time = (double)(PAP_stop-total_start)/1000.0; 

    removeFakes_TC();
    ocall_clock(&RF_stop); 
    RF_time = (double)(RF_stop - PAP_stop)/1000.0;
    //printf("Actual execution: PAP_time = %f, RF_time = %f\n", PAP_time, RF_time);
   
    // qsort is definitely faster!
    //mergeSort(outbuf, block_size, 0, n-1, &compare_keys);
    qsort(outbuf, N, block_size, compare);
    ocall_clock(&total_stop);
    ttime = (double)(total_stop-total_start)/1000.0;   
    return ttime;

  }

  #ifdef DETAILED_BOS_TIMING
    template <OSwap_Style oswap_style>
    double BucketOSort<oswap_style>::ORP(unsigned char *input, size_t N, size_t block_size, enc_ret *ret){
      long total_start, total_stop;
      double ttime;
      ocall_clock(&total_start);

      long PAP_stop, RF_stop;
      double PAP_time, RF_time;

      initializeBRN(); 
       
      processPacketsThroughBRN(input);

      flushBuffers();
      ocall_clock(&PAP_stop);
      PAP_time = (double)(PAP_stop-total_start)/1000.0; 

      removeFakes_TC();
      ocall_clock(&RF_stop); 
      RF_time = (double)(RF_stop - PAP_stop)/1000.0;
      //printf("Actual execution: PAP_time = %f, RF_time = %f\n", PAP_time, RF_time);
     
      ocall_clock(&total_stop);
      ttime = (double)(total_stop-total_start)/1000.0;   
      return ttime;
    }
  #else
    template <OSwap_Style oswap_style>
    double BucketOSort<oswap_style>::ORP(unsigned char *input, size_t N, size_t block_size, unsigned char *result_buf, enc_ret *ret) {
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Starting BORP (nthreads=%lu)\n", num_threads);
    #endif
      long total_start, total_stop, last_layer_start;
      double ttime;
      double lltime;
      ocall_clock(&total_start);

      //printf("Before initializeBRN()\n");
      initializeBRN(); 
       
      //printf("Before pptBRN()\n");
      threadpool_init(num_threads);

      //processPacketsThroughBRN(input);
      processPacketsThroughBRN_parallel(input, 0, N, num_threads);

      ocall_clock(&last_layer_start);
      ret->OSWAP_cb = OSWAP_COUNTER;
      flushBuffers();

      #ifndef BEFTS_MODE
        removeFakes_p1_parallel(0, b, num_threads);
      
        removeFakes_p2_parallel(0, b, result_buf, num_threads);
      #endif
      /*
      //printf("Before flushBuffers()\n");
      ocall_clock(&last_layer_start);
      flushBuffers();
    
      //printf("Before removeFakes_TC()\n");
      removeFakes_TC();
     
      ocall_clock(&total_stop);
      lltime = (double)(total_stop-last_layer_start)/1000.0;
      ret->last_layer_time=lltime;
      ttime = (double)(total_stop-total_start)/1000.0;   
      //printf("Done with ORP()\n");
      */
      threadpool_shutdown();
      ocall_clock(&total_stop);
      lltime = (double)(total_stop-last_layer_start)/1000.0;
      ret->last_layer_time=lltime;
      ttime = (double)(total_stop-total_start)/1000.0;
      ret->OSWAP_ap = OSWAP_COUNTER - ret->OSWAP_cb;
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Ending BORP(nthreads=%lu)\n", num_threads);
    #endif
      return ttime;
    }
  #endif


  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::calculatelog2f() {
    log2f = 0;
    int f_temp = 1;

    FOAV_SAFE2_CNTXT(calculatelog2f, f_temp, f)
    while(f_temp<f){
      f_temp=f_temp<<1;
      log2f+=1;
    }
  }

  template <OSwap_Style oswap_style>
  struct pptBRN_parallel_args {
      BucketOSort<oswap_style> *obj_ptr;
      unsigned char *input;
      size_t packet_start;
      size_t num_packets;
      size_t nthreads;
  };


  template<OSwap_Style oswap_style>
  void* BucketOSort<oswap_style>::processPacketsThroughBRN_parallel_launch(void *voidargs) {
      struct pptBRN_parallel_args<oswap_style> *args =
          (pptBRN_parallel_args<oswap_style> *) voidargs;
      (args->obj_ptr)->processPacketsThroughBRN_parallel(args->input, args->packet_start, args->num_packets, args->nthreads);
      return NULL;
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::processPacketsThroughBRN_parallel(unsigned char *input, 
  size_t packet_start, size_t num_packets, size_t nthreads) {

    FOAV_SAFE2_CNTXT(PPTBRN_parallel, packet_start, nthreads)
    if(nthreads<=1){
      // Execute the thread with processPacketsThroughBRN_MT() and the details
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting pptBRN_MT(N=%lu)\n", g_thread_id, num_packets);
    #endif
      processPacketsThroughBRN_MT(input, packet_start, num_packets);
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending pprBRN_MT(N=%lu)\n", g_thread_id, num_packets);
    #endif
    } else{ 
      // Part the problem over the nthreads this thread is responsible for
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting pptBRN_parallel(N=%lu, nthreads=%lu)\n", g_thread_id, num_packets, nthreads);
    #endif
      size_t lthreads = nthreads/2;
      size_t rthreads = nthreads - lthreads;
      size_t left_size = num_packets/2;
      size_t right_size = num_packets - left_size;

      unsigned char *input_left = input;
      unsigned char *input_right = input + (num_packets/2) * block_size;
      size_t packet_start_right = packet_start + left_size;

      threadid_t leftthreadid = g_thread_id + rthreads;
      BucketOSort<oswap_style> *obj_ptr = this;   
   
      pptBRN_parallel_args<oswap_style> leftargs = {
        obj_ptr, input, packet_start, left_size, lthreads
      };

      // dispatch leftthread with correct arguments.
      threadpool_dispatch(leftthreadid,
        processPacketsThroughBRN_parallel_launch,
        &leftargs);

      /* We will do the right half ourselves, on threads g_thread_id .. g_thread_id..rthreads-1 */
      processPacketsThroughBRN_parallel(input_right, packet_start_right, right_size, rthreads);
      threadpool_join(leftthreadid, NULL); 
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending pptBRN_parallel(N=%lu, nthreads=%lu)\n", g_thread_id, num_packets, nthreads);
    #endif
    }
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::processPacketsThroughBRN_MT(unsigned char *inbuf_ptr, 
        size_t packets_start, size_t num_packets) {
     
    size_t num_entry_msns = b/f;
    uint64_t packets_per_entry_msn = n/num_entry_msns;
    size_t msns_with_extra_packets = n % num_entry_msns; 
    size_t msn_start = packets_start / packets_per_entry_msn;
    size_t msn_stop = (packets_start + num_packets) / packets_per_entry_msn;
    size_t packets_left = num_packets;

    //printf("packets_per_entry_msn = %ld, num_packets = %ld\n", packets_per_entry_msn, num_packets);
    //printf("msn_start = %d, msn_stop = %d\n", msn_start, msn_stop);
    //printf("thread = %d, packets_start = %ld, num_packets = %ld\n", g_thread_id, packets_start, num_packets);
  
    const unsigned char *packet_ptr = inbuf_ptr;
    size_t ORP_packet_size = block_size + 16;
    unsigned char ORP_dummy_packet[ORP_packet_size];
    // Make this a dummy packet with dummy data
    generateORPDummy(ORP_dummy_packet);

    for(signed long i=msn_start; ((i<=msn_stop) && (packets_left!=0)); i++) {
      FOAV_SAFE2_CNTXT(PPTBRN_MT, i, msn_start)
      FOAV_SAFE2_CNTXT(PPTBRN_MT, i, msn_stop)
      FOAV_SAFE_CNTXT(PPTBRN_MT, packets_left)
      size_t packets_for_this_msn;
      

      if(i==msn_start) {
        FOAV_SAFE2_CNTXT(PPTBRN_MT, i, msn_start)
        size_t packets_msn_1 = packetsConsumedUptoMSN(i, msns_with_extra_packets, packets_per_entry_msn);
        size_t start_offset = packets_start - packets_msn_1;
        
        packets_for_this_msn = packets_per_entry_msn - start_offset + (i<msns_with_extra_packets);
        FOAV_SAFE2_CNTXT(PPTBRN_MT, packets_left, packets_for_this_msn)
        if(packets_left < packets_for_this_msn) {
          FOAV_SAFE2_CNTXT(PPTBRN_MT, packets_left, packets_for_this_msn)
          packets_for_this_msn = packets_left;
        }
        //printf("thread = %d, msn = %ld, packets_msn_1 = %ld, start_offset = %ld\n", g_thread_id, i, packets_msn_1, start_offset);
      }
      else{
        if(packets_left > (packets_per_entry_msn + (i<msns_with_extra_packets))) {
          FOAV_SAFE2_CNTXT(PPTBRN_MT, packets_left, packets_per_entry_msn)
          FOAV_SAFE2_CNTXT(PPTBRN_MT, i, msns_with_extra_packets)
          packets_for_this_msn = (packets_per_entry_msn + (i<msns_with_extra_packets));
        }
        else {
          packets_for_this_msn = packets_left;
        }
      }
      //printf("thread = %d, packets_per_entry_msn = %ld, packets_for_this_msn = %ld\n", g_thread_id, packets_per_entry_msn, packets_for_this_msn);
 
      //printf("msn = %ld of %ld, packets_for_this_msn = %ld\n", i, num_entry_msns, packets_for_this_msn); 
      for(unsigned int k=0; k<packets_for_this_msn; k++){
        FOAV_SAFE2_CNTXT(PPTBRN_MT, k, packets_for_this_msn)
        // New ORP Packet to be generated
        unsigned char new_ORP_packet[ORP_packet_size];
        // Reset dummy packet (in case it got overwritten with other packets in the previous 
        // processAndPropogate pass
        setORPDummy(ORP_dummy_packet);
        unsigned char *new_packet_ptr = new_ORP_packet;

        // Give the packet a random label for ORP
        uint16_t random_label;
    
        getRandomBytes((unsigned char*) &random_label, sizeof(uint16_t));
        random_label = random_label % b;
        //bucket_distribution[random_label]++;
        ((uint64_t*) new_packet_ptr)[1] = random_label;
        // No need to set eviction_stream for ORP packet, that will be assigned by the MSN 
        // that processes this packet 
        // Shift new_packet_ptr to account for the 2 added fields for ORP packets
        new_packet_ptr = new_packet_ptr+(2*sizeof(uint64_t));
        memcpy(new_packet_ptr, packet_ptr, block_size); 
        //printf("Processing packet %ld\n", *((uint64_t*) new_packet_ptr));
        
        //printf("ProcessPackets invoked with real packet label = %d\n", random_label);
        nodes[0][i]->processAndPropogate(new_ORP_packet);

        // Process a dummy packet
        nodes[0][i]->processAndPropogate(ORP_dummy_packet);
        packet_ptr+=block_size;
      }
      packets_left-=packets_for_this_msn;
    }
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::flushBuffers_MT(int depth, int msn_start, int msn_stop)  {
    for(int j=msn_start; j<msn_stop; j++) {
      for(unsigned int k=0; k< B*f; k++)  {
        FOAV_SAFE2_CNTXT(flushBuffers_MT, B, f)
        FOAV_SAFE2_CNTXT(flushBuffers_MT, j, k)
        FOAV_SAFE2_CNTXT(flushBuffers_MT, msn_start, msn_stop)
        unsigned char dummy_packet[ORP_packet_size];
        setORPDummy(dummy_packet);
        nodes[depth][j]->processAndPropogate(dummy_packet);
      }
    }
  }

  template <OSwap_Style oswap_style>
  struct removeFakes_p1_parallel_args {
      BucketOSort<oswap_style> *obj_ptr;
      int bucket_start;
      int bucket_stop;
      size_t nthreads;
  };

  template<OSwap_Style oswap_style>
  void* BucketOSort<oswap_style>::removeFakes_p1_parallel_launch(void *voidargs) {
      struct removeFakes_p1_parallel_args<oswap_style> *args =
          (removeFakes_p1_parallel_args<oswap_style> *) voidargs;
      (args->obj_ptr)->removeFakes_p1_parallel(args->bucket_start, args->bucket_stop, args->nthreads);
      return NULL;
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes_p1_parallel(int bucket_start, int bucket_stop, size_t nthreads) {
    // If we are at a single thread, or if there is only one bucket we stop and let that thread execute removeFakes_MT_p1
    FOAV_SAFE2_CNTXT(removeFakes_p1, nthreads, bucket_start)
    FOAV_SAFE_CNTXT(removeFakes_p1, bucket_stop)
    if(nthreads<=1 || (bucket_stop == bucket_start)) {
      FOAV_SAFE2_CNTXT(removeFakes_p1, nthreads, bucket_start)
      FOAV_SAFE_CNTXT(removeFakes_p1, bucket_stop)
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting removeFakes_MT_p1(buckets=%lu)\n", g_thread_id, 1);
    #endif
      removeFakes_p1(bucket_start, bucket_stop);
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending removeFakes_MT_p1(buckets=%lu)\n", g_thread_id, 1);
    #endif
    }
    else{
      // Thread dispatch
      // Part the problem over the nthreads this thread is responsible for
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting removeFakes_MT_p1(buckets=%lu, nthreads=%lu)\n", g_thread_id, (bucket_stop-bucket_start), nthreads);
    #endif
      size_t lthreads = nthreads/2;
      size_t rthreads = nthreads - lthreads;
      int mid_bucket = (bucket_stop + bucket_start)/2;
      int left_start = bucket_start;
      int left_stop = mid_bucket;
      int right_start = mid_bucket;
      int right_stop = bucket_stop;

      threadid_t leftthreadid = g_thread_id + rthreads;
      BucketOSort<oswap_style> *obj_ptr = this;   
   
      removeFakes_p1_parallel_args<oswap_style> leftargs = {
        obj_ptr, left_start, left_stop, lthreads
      };

      // dispatch leftthread with correct arguments.
      threadpool_dispatch(leftthreadid,
        removeFakes_p1_parallel_launch,
        &leftargs);

      /* We will do the right half ourselves, on threads g_thread_id .. g_thread_id..rthreads-1 */
      removeFakes_p1_parallel(right_start, right_stop, rthreads);
      threadpool_join(leftthreadid, NULL); 
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending removeFakes_MT_p1_parallel(buckets=%lu, nthreads=%lu)\n", g_thread_id, (bucket_stop-bucket_start), nthreads);
    #endif
    }

  }

  /*
    Go through each bucket of outbuf and populate the count of real packets in each of them
    to the BOS's local array reals_in_bucket
  */
  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes_p1(int bucket_start, int bucket_stop) {
    // Calculate outbuf offset based on bucket start
    // Iterate over buckets bucket_start to bucket_stop and update the count of 
    // real_packets to reals_in_bucket

    for(int i = bucket_start; i < bucket_stop; i++) {
      FOAV_SAFE2_CNTXT(removeFakes_p1, i, bucket_start)
      FOAV_SAFE_CNTXT(removeFakes_p1, bucket_stop)
      unsigned char *bkt_ptr = outbuf + i * (uint64_t)(num_packets_in_outbuf_bucket * block_size);
      bool *selected_ptr = selected_in_bucket + (i*num_packets_in_outbuf_bucket);
      uint32_t num_reals = 0;
      // Scan over the contents of the bucket and track the number of reals in it obliviously    
      for(size_t k=0; k<num_packets_in_outbuf_bucket; k++) {
        FOAV_SAFE2_CNTXT(removeFakes_p1, k, num_packets_in_outbuf_bucket)
        uint32_t real_flag = !(isDummy(bkt_ptr));
        num_reals+=real_flag;
        bkt_ptr+=block_size;
        selected_ptr[k]=real_flag;
      }
      //printf("removeFakes_p1: num_reals for bucket %d = %d\n", i, num_reals);
      reals_in_bucket[i] = num_reals;
    }
  }

  template <OSwap_Style oswap_style>
  struct removeFakes_p2_parallel_args {
      BucketOSort<oswap_style> *obj_ptr;
      int bucket_start;
      int bucket_stop;
      unsigned char *result_buf;
      size_t nthreads;
  };

  template<OSwap_Style oswap_style>
  void* BucketOSort<oswap_style>::removeFakes_p2_parallel_launch(void *voidargs) {
      struct removeFakes_p2_parallel_args<oswap_style> *args =
          (removeFakes_p2_parallel_args<oswap_style> *) voidargs;
      (args->obj_ptr)->removeFakes_p2_parallel(args->bucket_start, args->bucket_stop,
         args->result_buf, args->nthreads);
      return NULL;
  }

  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes_p2_parallel(int bucket_start, int bucket_stop, unsigned char *result_buf, size_t nthreads) {
    // If we are at a single thread, or if there is only one bucket we stop and let that thread execute removeFakes_MT_p1
    FOAV_SAFE_CNTXT(removeFakes_p2, bucket_stop)
    FOAV_SAFE2_CNTXT(removeFakes_p2, nthreads, bucket_start)
    if(nthreads<=1 || (bucket_stop == bucket_start)) {
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting removeFakes_MT_p2(buckets=%lu)\n", g_thread_id, nthreads);
    #endif
      removeFakes_p2(bucket_start, bucket_stop, result_buf, nthreads);
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending removeFakes_MT_p2(buckets=%lu)\n", g_thread_id, nthreads);
    #endif
    }
    else{
      // Thread dispatch
      // Part the problem over the nthreads this thread is responsible for
    #ifdef VERBOSE_TIMINGS_BORP
      unsigned long start = printf_with_rtclock("Thread %u starting removeFakes_MT_p2(buckets=%lu, nthreads=%lu)\n", g_thread_id, (bucket_stop-bucket_start), nthreads);
    #endif
      size_t lthreads = nthreads/2;
      size_t rthreads = nthreads - lthreads;
      int mid_bucket = (bucket_stop + bucket_start)/2;
      int left_start = bucket_start;
      int right_stop = bucket_stop;

      threadid_t leftthreadid = g_thread_id + rthreads;
      BucketOSort<oswap_style> *obj_ptr = this;   
   
      removeFakes_p2_parallel_args<oswap_style> leftargs = {
        obj_ptr, left_start, mid_bucket, result_buf, lthreads
      };

      // dispatch leftthread with correct arguments.
      threadpool_dispatch(leftthreadid,
        removeFakes_p2_parallel_launch,
        &leftargs);

      /* We will do the right half ourselves, on threads g_thread_id .. g_thread_id..rthreads-1 */
      removeFakes_p2_parallel(mid_bucket, right_stop, result_buf, rthreads);
      threadpool_join(leftthreadid, NULL); 
    #ifdef VERBOSE_TIMINGS_BORP
      printf_with_rtclock_diff(start, "Thread %u ending removeFakes_MT_p2_parallel(buckets=%lu, nthreads=%lu)\n", g_thread_id, (bucket_stop-bucket_start), nthreads);
    #endif
    }
  }


  template <OSwap_Style oswap_style>
  void BucketOSort<oswap_style>::removeFakes_p2(int bucket_start, int bucket_stop, unsigned char *result_buf, size_t nthreads) {
    //printf("In removeFakes_MT_p2, thread = %d, outbuf = %p, bucket_start = %d, bucket_stop = %d\n", g_thread_id, outbuf, bucket_start, bucket_stop); 
    FOAV_SAFE2_CNTXT(removeFakes_p2, bucket_stop, bucket_start)
    uint64_t i =0;
    FOAV_SAFE_CNTXT(removeFakes_p2, i)
    size_t offset =  (num_packets_in_outbuf_bucket * bucket_start * block_size);
    unsigned char *outbuf_ptr = outbuf + offset;
    //printf("outbuf = %p, num_packets_in_outbuf_bucket = %u, bucket_start = %u, block_size = %u, offset = %u, outbuf_ptr = %p\n", outbuf, num_packets_in_outbuf_bucket, bucket_start, block_size, offset, outbuf_ptr);
    size_t result_buf_offset=0;
    FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_start)
    for(i =0; i<bucket_start; i++){
      FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_start)
      result_buf_offset+=(reals_in_bucket[i]);
    }
    FOAV_SAFE2_CNTXT(removeFakes_p2, bucket_stop, bucket_start)
    result_buf_offset*=block_size;
    unsigned char *result_ptr = result_buf+result_buf_offset;
    //printf("result_buf = %x, result_buf_offset = %ld, result_ptr = %x\n", result_buf, result_buf_offset, result_ptr);
     
    FOAV_SAFE2_CNTXT(removeFakes_p2, result_ptr, result_buf_offset)
    FOAV_SAFE2_CNTXT(removeFakes_p2, bucket_stop, bucket_start)
    FOAV_SAFE2_CNTXT(removeFakes_p2, result_buf, nthreads)
    FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_start)
    for(i=bucket_start; i<bucket_stop; i++){
      FOAV_SAFE2_CNTXT(removeFakes_p2, result_buf, nthreads)
      FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_start)
      FOAV_SAFE_CNTXT(removeFakes_p2, bucket_stop)
      size_t reals_in_this_bucket = reals_in_bucket[i];
      //printf("in RFMT_p2, reals_in_this_bucket = %ld\n", reals_in_this_bucket);
      bool *selected_ptr = selected_in_bucket + (i*num_packets_in_outbuf_bucket);

      // Invoke TightCompaction on each bucket
      TightCompact_v2(outbuf_ptr, num_packets_in_outbuf_bucket, block_size, selected_ptr);

      // NOTE: Reusing the portion of selected_list from TC call, since that is no longer useful, and MarkHalf should overwrite it correctly.
      //       We also allow RS_M2 to parallelize further if the thread still owns subthreads. (We don't do this in removeFakes_p1)
      // We reset <oswap_style> here since <oswap_style> for BORP is chosen based on the ORP_packet_size = block_size + 16 
      FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
      if(block_size==4){
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
        RecursiveShuffle_M2_inner_parallel<OSWAP_4>(outbuf_ptr, reals_in_this_bucket, block_size, selected_ptr, nthreads);
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
      } else if(block_size==8){
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
        RecursiveShuffle_M2_inner_parallel<OSWAP_8>(outbuf_ptr, reals_in_this_bucket, block_size, selected_ptr, nthreads);
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
      } else if(block_size%16==0) {
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
        RecursiveShuffle_M2_inner_parallel<OSWAP_16X>(outbuf_ptr, reals_in_this_bucket, block_size, selected_ptr, nthreads);
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
      } else {
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
        RecursiveShuffle_M2_inner_parallel<OSWAP_8_16X>(outbuf_ptr, reals_in_this_bucket, block_size, selected_ptr, nthreads);
        FOAV_SAFE_CNTXT(removeFakes_p2, block_size)
      }
 
      memmove(result_ptr, outbuf_ptr, reals_in_this_bucket * block_size);
      outbuf_ptr = outbuf_ptr+ (block_size * num_packets_in_outbuf_bucket);
      result_ptr += (block_size * reals_in_this_bucket);
      FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_start)
      FOAV_SAFE2_CNTXT(removeFakes_p2, result_buf, nthreads)
      FOAV_SAFE2_CNTXT(removeFakes_p2, i, bucket_stop)
      FOAV_SAFE2_CNTXT(removeFakes_p2, result_ptr, result_buf_offset)
      FOAV_SAFE_CNTXT(removeFakes_p2, outbuf_ptr)
    }  
  }
#endif
