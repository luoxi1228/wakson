
#include "MergeSplitNode.hpp"

int process_MSN(unsigned char *packet, bos_params *params, size_t repeat) {
  size_t block_size = params->ORP_packet_size;

  //ORP_packet_size will never be 8!
  if(block_size%16==0){
    FOAV_SAFE_CNTXT(proces_MSN, block_size)
    MergeSplitNode<OSWAP_16X> *tuner_MSN = new MergeSplitNode<OSWAP_16X>(params, false);
    for(size_t i = 0; i<repeat; i++){
      FOAV_SAFE2_CNTXT(process_MSN, i, repeat)
      tuner_MSN->process(packet);
    }
    delete tuner_MSN;    
  } /*else {
    MergeSplitNode<OSWAP_8_16X> *tuner_MSN = new MergeSplitNode<OSWAP_8_16X>(params);
    for(size_t i = 0; i<repeat; i++){
      tuner_MSN->process(packet);
    }
    delete tuner_MSN;
  }*/
  

  return 0;
}
