#ifndef __PROCESSING_NODE__
  #define __PROCESSING_NODE__

  #include "Globals_Enclave.hpp"
  #include "../Globals.hpp"
  #include "../CONFIG.h"
  #include "oasm_lib.h"

  class ProcessingNode{
    public:
      unsigned int f;
      unsigned int B;
      unsigned int occupancy;
      
      // the MSN needs to know how many layers, and which layer of the BRN it is in so that it
      // can correctly route by looking at right log_2f bits of the destination
      // node_in_layer has no utility for MSN, it doesn't need to know which node in
      // a given layer it is, but this could be useful to debug processAndPropogate()
      unsigned int d;
      unsigned int layer;
      unsigned int node_in_layer;

      unsigned int log2f;
      unsigned int mask;

      unsigned int evict_start_param;
      // current_evict_stream will round robin over 1 to f  
      uint64_t current_evict_stream;    

      // When used in Bucket Oblivious Sort, this is the handle to the 
      // ProcessingNodes this node is connected to in the next layer
      ProcessingNode **eviction_streams; 

      // ProcessingNode Inherited Virtual Functions/Interfaces:      
      virtual void initializeBuffer()=0; 
      virtual int setEvictionStreams(ProcessingNode **nodes)=0;
      virtual ProcessingNode**  getEvictionStreams()=0;
      virtual void setOutbufPtr(unsigned char *outbuf_ptr)=0;

      virtual int getBufferOccupancy()=0;
      /* Process a packet and then call processAndPropogate for the next node (based
          on the current_evict_stream).
      */
      virtual void processAndPropogate(unsigned char* packet_in_out)=0;

      // Debug/display functions
      virtual void printNodePosition()=0;
  }; 

#endif
