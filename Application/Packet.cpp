
#include "Packet.hpp"

void fillPacketData(packet* packet_in){
  packet_in->data = (unsigned char*) malloc(DATA_SIZE);
  char curr_char = 'A';
  for(int i=0; i<DATA_SIZE;i++){
    (packet_in->data)[i]=curr_char;
    curr_char='A'+((i+1)%26);
  }
}

packet* generate_packet(uint64_t n, uint64_t b){
  packet *new_packet = (packet *) malloc(sizeof(packet));
  uint64_t r = rand()%b;
  //printf("%ld\n", r);

  new_packet->in_address=0;
  fillPacketData(new_packet);
  new_packet->out_address=r;

  return new_packet;
}

packet* generate_packet(unsigned int l, unsigned int k, unsigned int f) {
  packet *new_packet = (packet *) malloc(sizeof(packet));
  unsigned int r = rand()%k;

  uint64_t destination=UINT64_MAX;
  //TODO: Fix in_address, not bothered with it for now.
  new_packet->in_address=0;

  if(r<l) {
    //generate real packet
    // +1 to shift it to 1 to f as the stream indices
    // 
    // When packet is being submitted:
    // out_address = -1 for dummy
    // 
    // When enclave returns a packet:
    // out_address =  0 will be used by the enclave to signal no_packets evicted
    // out_address = -1 for Overflow
    destination = (rand()%f)+1; 
  }
  else {
    //generate dummy packet
  }  
  fillPacketData(new_packet);

  new_packet->out_address=destination;
  return new_packet;
}

unsigned char *serialize_packet(packet *packet_in) {
  unsigned char *serialized_packet = (unsigned char*) malloc(SERIALIZED_PACKET_SIZE); 
  unsigned char *packet_ptr = serialized_packet;

  memcpy(packet_ptr, (unsigned char*) &(packet_in->in_address), sizeof(uint64_t));
  packet_ptr+=sizeof(uint64_t);
  
  memcpy(packet_ptr, (unsigned char*) &(packet_in->out_address), sizeof(uint64_t));
  packet_ptr+=sizeof(uint64_t);

  if(packet_in->data==NULL)
    printf("Uhoh, tried to serialize a packet without data.\n");
  else
    memcpy(packet_ptr, packet_in->data, DATA_SIZE);
  
  return serialized_packet;
}

packet* deserialize_packet(unsigned char *serialized_packet) {
  packet *new_packet = (packet*) malloc(sizeof(packet));
  new_packet->data = (unsigned char*) malloc(DATA_SIZE);

  unsigned char *packet_ptr = serialized_packet;
 
  memcpy((unsigned char*) &(new_packet->in_address), packet_ptr, sizeof(uint64_t));
  packet_ptr+=sizeof(uint64_t);

  memcpy((unsigned char*) &(new_packet->out_address), packet_ptr, sizeof(uint64_t));
  packet_ptr+=sizeof(uint64_t);

  memcpy(new_packet->data, packet_ptr, DATA_SIZE);

  return new_packet;
}


//NOTE: Convenience hijacking the last byte of data with '\0' 
void displayPacket(packet *packet_in) {
  char data[DATA_SIZE+1];
  memcpy(data, packet_in->data, DATA_SIZE);  
  data[DATA_SIZE]='\0';
  printf("(In = %ld, Out = %ld) : %s\n", packet_in->in_address, packet_in->out_address, data);
}

void free_packet(packet *packet_in) {
  free(packet_in->data);
  free(packet_in);
}
