#include "Bucket.h"
#include <stdexcept>

Bucket::Bucket() : subtree_min(BlockHeader()) {
    for (int i=0; i < BUCKET_CAPACITY; ++i) {
        this->blocks[i] = Block();
    }
}

void Bucket::addBlock(const Block& new_blk) {
  for (int i = 0; i < BUCKET_CAPACITY; i++) {
    if (this->blocks[i].h.t == -1) {
        this->blocks[i] = new_blk;
        break;
    }
  }
}

std::vector<Block> Bucket::removeBlock(long t) {
  std::vector<Block> ret;
  for (int i = 0; i < BUCKET_CAPACITY; i++) {
    if (this->blocks[i].h.t == t) {
      ret.push_back(this->blocks[i]);
      this->blocks[i].h.t = -1;
      break;
    }
  }
  return ret;
}
