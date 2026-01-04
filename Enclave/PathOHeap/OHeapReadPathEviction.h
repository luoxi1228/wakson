#ifndef POHEAP_OHEAPREADPATHEVICTION_H
#define POHEAP_OHEAPREADPATHEVICTION_H

#include "OHeapInterface.h"
#include "RandForOHeapInterface.h"
#include "UntrustedStorageInterface.h"

class OHeapReadPathEviction : public OHeapInterface {
 public:
  OHeapReadPathEviction(UntrustedStorageInterface* storage,
      RandForOHeapInterface* rand_gen, int num_blocks, bool type_hiding_security);

  BlockHeader findMin();
  
  /** Returns the leaf_id (pos) and t of the inserted block. */
  std::pair<int, int> insertBlock(int k, const unsigned char *v);
  
  /** Returns the deleted block if one matching the parameters exists, otherwise a dummy block. */
  Block deleteBlock(int leaf, long t);
  
  Block extractMin();

  /** Ill-defined if no block with matching leaf and t exists. */
  std::pair<int, int> decreaseKey(int leaf, long t);
  std::pair<int, int> increaseKey(int leaf, long t);

  std::vector<Block>& getStash();
  int getStashSize();
  int getNumLeaves();
  int getNumLevels();
  int getNumBlocks();
  int getNumBuckets();
  // friend std::ostream& operator<<(std::ostream& out, const OHeapReadPathEviction& oheap);

 private:
  int num_levels;
  int num_leaves;
  int num_blocks;
  int num_buckets;
  long t;  // timestamp or number of ops performed so far
  bool type_hiding_security;

  UntrustedStorageInterface* storage;
  RandForOHeapInterface* rand_gen;
  // The stash will be allocated just once, at object constructor time.  The
  // stash is usually quite small (~5 elements for a million-item tree), but
  // during eviction, all the blocks on a path are temporarily added to it
  // before being reinserted back into the path.  So we set two stash sizes:
  // the short_stashsize, which is the size outside of eviction, and the
  // long_stashsize, which is the size during eviction.  We allocate stash to
  // be the long_stashsize, but we only use the first short_stashsize elements
  // except during eviction.  When eviction is done, we will move all of the
  // non-dummy elements (obliviously) to the first short_stashsize positions.
  std::vector<Block> stash;
  size_t short_stashsize;
  size_t long_stashsize;

  BlockHeader onlyFindMin();
  Block onlyDeleteBlock(int leaf, long t);
  std::pair<int, int> onlyInsertBlock(int k, const unsigned char *v);
  void insertDummyBlock();
  void evictAfterInsert();

  /**
   * Evicts the path corresponding to the given leaf, then runs updateMin.
   * (It never makes sense to evict without updating min.)
   */
  void evictAndUpdateMin(int leaf);
  void evictAndUpdateMin(int leaf, std::vector<Bucket>& path, const std::vector<Bucket>& pathsiblings);

  /**
   * Each leaf corresponds to a specific path made up of one bucket per level; in reality buckets
   * are placed linearly in storage. Returns the correct storage index so that bucket can be
   * retrived.
   */
  int P(int leaf, int level);

  /**
   * Operation P.ReadNRm from the OHeap paper where leaf identifies the path.
   * Returns the removed Block if one matching the parameters exists, otherwise a dummy block.
   */
  Block readAndRemove(int leaf, long t);

  /** Operation P.updateMin from the OHeap paper where leaf identifies the path. */
  void updateMin(int leaf, std::vector<Bucket>& path, const std::vector<Bucket>& pathsiblings);
};

#endif 
