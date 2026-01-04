#ifndef POHEAP_SERVERSTORAGE_H
#define POHEAP_SERVERSTORAGE_H

#include "UntrustedStorageInterface.h"

class ServerStorage : public UntrustedStorageInterface {
 public:
  ServerStorage();

  /** Must be called before reading or writing any buckets. Cannot be changed after being set. */
  void setCapacity(int total_num_of_buckets);

  /** Returns a reference to the root of the tree. */
  Bucket& readRoot();

  /** Returns a copy of the read bucket. The original bucket cannot be modified directly. */
  Bucket readBucket(int position);

  /** Copies the given bucket into this storage. */
  void writeBucket(int position, const Bucket& bucket_to_write);

  /** Reads the path to the given leaf into the local path, as well as the
   * siblings of all the (non-root) elements on that path, and returns a
   * reference to the path, and the vector of siblings. */
  std::pair<std::vector<Bucket>&,const std::vector<Bucket>&> readPath(int leaf);

  /** Write back the local path into storage, invalidating the local path. */
  void writebackPath();

private: 
  static bool is_initialized;
  static bool is_capacity_set;
  int capacity;
  int num_levels;
  int num_leaves;
  std::vector<Bucket> buckets;
  std::vector<Bucket> curpath;
  std::vector<Bucket> pathsiblings;
  int curleaf;

  void checkRWValid(int position);
};

#endif
