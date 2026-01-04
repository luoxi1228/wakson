
#ifndef __PACKET__

  #include <cstdio>
  #include <string.h>
  #include <stdlib.h>
  #include <stdint.h>
  #include <random>
  #include "../CONFIG.h"
  #include "../Globals.hpp"

  packet* generate_packet(unsigned int l, unsigned int k, unsigned int f);
  packet* generate_packet(uint64_t n, uint64_t b);
  unsigned char *serialize_packet(packet *packet_in);
  packet* deserialize_packet(unsigned char *serialized_packet);
  void fillPacketData(packet *packet_in);
  void displayPacket(packet *packet_in);
  void free_packet(packet *packet_in);

  #define __PACKET__
#endif
