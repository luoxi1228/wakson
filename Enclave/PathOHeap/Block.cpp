#include "Block.h"
#include <climits>

BlockHeader::BlockHeader(int leaf_id, long t, int k) : leaf_id(leaf_id), t(t), k(k) {}

BlockHeader::BlockHeader() : leaf_id(-1), t(-1), k(INT_MAX) {}

Block::Block(int leaf_id, long t, int k, const unsigned char *v) : h(leaf_id, t, k) {
  memmove(this->v, v, BLOCK_SIZE);
}

Block::Block() {
  memset(this->v, 0, BLOCK_SIZE);
}

bool operator<(const Block& b1, const Block& b2) {
  return (b1.h.k < b2.h.k) || (b1.h.k == b2.h.k && b1.h.t < b2.h.t);
}

bool operator==(const Block& b1, const Block& b2) {
  return !memcmp(&b1.h, &b2.h, sizeof(BlockHeader));
}

bool operator!=(const Block& b1, const Block& b2) {
  return !(b1 == b2);
}

/*
std::ostream& operator<<(std::ostream& out, const Block& b) {
  if (b.leaf_id != -1) {
    out << "Block (leaf, t, k, v) = " << b.leaf_id << " " << b.t << " " << b.k << " [ ";
    for (std::byte byte : b.v) {
      out << std::to_integer<int>(byte) << " ";
    }
    out << "]";
  } else {
    out << "Dummy block";
  }
  return out;
}
*/
