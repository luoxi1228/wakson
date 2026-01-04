#ifndef POHEAP_BLOCK_H
#define POHEAP_BLOCK_H

#include <array>
#include <cstddef>
#include <iostream>

struct BlockHeader {
  int t;  // tau or timestamp
  int k;
  int leaf_id;  // aka pos
  unsigned int heapno;

  BlockHeader(int leaf_id, long t, int k);
  BlockHeader();
};

class Block {
 public:
  static const int BLOCK_SIZE = 288;

  BlockHeader h;
  unsigned char v[BLOCK_SIZE];

  Block(int leaf_id, long t, int k, const unsigned char *v);
  Block();  /** dummy index */
};

bool operator<(const Block& b1, const Block& b2);
bool operator==(const Block& b1, const Block& b2);
bool operator!=(const Block& b1, const Block& b2);
//std::ostream& operator<<(std::ostream& out, const Block& b);

#endif
