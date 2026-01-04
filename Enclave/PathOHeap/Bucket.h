#ifndef POHEAP_BUCKET_H
#define POHEAP_BUCKET_H

#include "Block.h"
#include <vector>

class Bucket {
 public:
  static const int BUCKET_CAPACITY = 3;

  BlockHeader subtree_min;
  Block blocks[BUCKET_CAPACITY];

  Bucket();

  /**
   * Caller is responsible for making sure new_blk is allowed in this bucket!
   * Caller is responsible for updating this bucket's subtree_min afterward!
   * Copies the given block into this bucket.
   */
  void addBlock(const Block& new_blk);

  /**
   * Caller is responsible for updating this bucket's subtree_min afterward!
   * Returns the removed block if a timestamp match is found.
   */
  std::vector<Block> removeBlock(long t);

};

#endif
