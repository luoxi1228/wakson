#ifndef POHEAP_SGXRAND_H
#define POHEAP_SGXRAND_H

#include "RandForOHeapInterface.h"

class SGXRand : public RandForOHeapInterface {
 public:
  SGXRand();

  /** Returns a random integer between 0 (inclusive) and bound (exclusive) */
  int getRandomLeaf();

  void setBound(int num_leaves);

 private:
  unsigned int maxval;
  unsigned int mask;
};

#endif
