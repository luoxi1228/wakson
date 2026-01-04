#include "OHeapReadPathEviction.h"
#include <cmath>
#include <algorithm>
#include <climits>
#include "utils.hpp"

OHeapReadPathEviction::OHeapReadPathEviction(UntrustedStorageInterface* storage,
    RandForOHeapInterface* rand_gen, int num_blocks, bool type_hiding_security)
    : storage(storage), rand_gen(rand_gen), num_blocks(num_blocks),
      type_hiding_security(type_hiding_security) {
  printf("sizeof(Block) = %d\nsizeof(Bucket) = %d\n", sizeof(Block), sizeof(Bucket));
  this->num_levels = ceil(log2(num_blocks)) + 1;
  this->num_buckets = (1<<num_levels) - 1;
  if (this->num_buckets*Bucket::BUCKET_CAPACITY < this->num_blocks) {  // potential precision loss
    throw std::runtime_error("Not enough space for the acutal number of blocks.");
  }
  this->num_leaves = 1<<(num_levels - 1);
  this->t = 0;

/*
  std::cout << "Initializing OHeap with " << this->num_levels << " levels, " <<
      this->num_blocks << " blocks, " << this->num_buckets << " buckets" << std::endl;
*/

  this->rand_gen->setBound(num_leaves);
  this->storage->setCapacity(num_buckets);

  // These constants may use some tweaking.  If they're too low, you may get an
  // exception thrown if the stash fills (which will appear as an illegal
  // instruction error from the enclave if uncaught).  If too high, all
  // operations on the stash will take (linearly) longer.
  this->short_stashsize = (this->num_levels)/3;
  this->long_stashsize = 2*(this->num_levels);
  this->stash.assign(this->long_stashsize, Block());
}

BlockHeader OHeapReadPathEviction::findMin() {
  BlockHeader ret = this->onlyFindMin();

  if (this->type_hiding_security) {
    this->onlyDeleteBlock(-1, -1);
    this->insertDummyBlock();
  }

  return ret;
}

Block OHeapReadPathEviction::deleteBlock(int leaf, long t) {
  if (this->type_hiding_security) {
    this->onlyFindMin();
  }

  Block ret = this->onlyDeleteBlock(leaf, t);

  if (this->type_hiding_security) {
    this->insertDummyBlock();
  }

  this->t++;
  return ret;
}

std::pair<int, int> OHeapReadPathEviction::insertBlock(int k, const unsigned char *v) {
  if (this->type_hiding_security) {
    this->onlyFindMin();
    this->onlyDeleteBlock(-1, -1);
  }

  auto ret = std::pair<int, int>(-1, -1);

  if (k != INT_MAX) {
    ret = this->onlyInsertBlock(k, v);
  } else {
    this->insertDummyBlock();
  }

  this->t++;
  return ret;
}

Block OHeapReadPathEviction::extractMin() {
  BlockHeader min_block = this->onlyFindMin();
  Block ret = this->onlyDeleteBlock(min_block.leaf_id, min_block.t);

  if (this->type_hiding_security) {
    this->insertDummyBlock();
  }

  this->t++;
  return ret;
}

std::pair<int, int> OHeapReadPathEviction::decreaseKey(int leaf, long t) {
  if (this->type_hiding_security) {
    this->onlyFindMin();
  }

  Block old_block = this->onlyDeleteBlock(leaf, t);
  auto ret = this->onlyInsertBlock(old_block.h.k - 1, old_block.v);

  this->t++;
  return ret;
}

std::pair<int, int> OHeapReadPathEviction::increaseKey(int leaf, long t) {
  if (this->type_hiding_security) {
    this->onlyFindMin();
  }

  Block old_block = this->onlyDeleteBlock(leaf, t);
  auto ret = this->onlyInsertBlock(old_block.h.k + 1, old_block.v);

  this->t++;
  return ret;
}

// An oblivious constant-time implementation of inserting the given block
// into the stash.  It is allowed for the given block to be a dummy (in
// which case this is a no-op).  Returns 0 on success (including the
// inserting a dummy case), -1 if the stash is full.
static int stashInsert(Block *stash, size_t stashsize, const Block &block)
{
    // If we're given a dummy block, initialize ret to 0, otherwise -1
    int ret = -(block.h.t != -1);
    Block *stashend = stash+stashsize;
    while (stash < stashend) {
        const Block *blockptr = &block;
        __asm__ (
            "# If %[flag] is 1, move 304 bytes (38*8) from %[src] to %[dst] and set %[ret] to 0\n"
            "movl $7, %%ecx\n"          // Loop 7 times moving 40 bytes at a time = 280 bytes
            ".stashInsertL%=:\n"
            "test %[flag], %[flag]\n"
            "movq (%[dst]), %%rax\n"
            "movq 8(%[dst]), %%rbx\n"
            "movq 16(%[dst]), %%rdi\n"
            "movq 24(%[dst]), %%r8\n"
            "movq 32(%[dst]), %%r9\n"
            "cmovneq (%[src]), %%rax\n"
            "cmovneq 8(%[src]), %%rbx\n"
            "cmovneq 16(%[src]), %%rdi\n"
            "cmovneq 24(%[src]), %%r8\n"
            "cmovneq 32(%[src]), %%r9\n"
            "movq %%rax, (%[dst])\n"
            "movq %%rbx, 8(%[dst])\n"
            "movq %%rdi, 16(%[dst])\n"
            "movq %%r8, 24(%[dst])\n"
            "movq %%r9, 32(%[dst])\n"
            "addq $40, %[dst]\n"
            "addq $40, %[src]\n"
            "dec %%ecx\n"
            "jne .stashInsertL%=\n"
            "test %[flag], %[flag]\n"
            "cmovnel %%ecx, %[ret]\n"   // Move 0 into %[ret]
            "movq (%[dst]), %%rax\n"    // Move the last 24 bytes
            "movq 8(%[dst]), %%rbx\n"
            "movq 16(%[dst]), %%r8\n"
            "cmovneq (%[src]), %%rax\n"
            "cmovneq 8(%[src]), %%rbx\n"
            "cmovneq 16(%[src]), %%r8\n"
            "movq %%rax, (%[dst])\n"
            "movq %%rbx, 8(%[dst])\n"
            "movq %%r8, 16(%[dst])\n"
            "addq $24, %[dst]\n"
            : [ret] "+r" (ret), [dst] "+r" (stash), [src] "+r" (blockptr)
            : [flag] "r" ((stash->h.t == -1) & ret)
            : "cc", "memory", "rax", "rbx", "rcx", "rdi", "r8", "r9"
        );
    }
    return ret;
}

// An oblivious constant-time implementation of extracting a block from the
// stash.  Moves any one block in the stash whose leaf_id matches the given
// one up to the given mask into dst, and marks the moved block as a dummy
// in the stash.  If no block in the stash matches the leaf_id (up to the
// mask), a dummy block is written to dst.
static void stashExtract(Block *stash, size_t stashsize,
        int leaf_id, unsigned int mask, Block *dst)
{
    dst->h.t = -1;
    int seen = 0;  // seen is set to 1 if a block is extracted, and no other
                   // block will be extracted once seen is 1
    Block *stashend = stash+stashsize;
    while (stash < stashend) {
        Block *dstblock = dst;
        __asm__ (
            "# If %[flag] is 0, move 304 bytes (38*8) from %[src] to %[dst]\n"
            "# and set %[seen] to 1 and *%[src] (== stash->h.t) to -1\n"
            "xor %%edi, %%edi\n"
            "dec %%edi\n"               // edi <- -1
            "test %[flag], %[flag]\n"
            "movq (%[dst]), %%rax\n"    // Move the first 24 bytes and set *%[src] to -1
            "movq 8(%[dst]), %%rbx\n"
            "movq 16(%[dst]), %%r8\n"
            "cmoveq (%[src]), %%rax\n"
            "movl (%[src]), %%r9d\n"
            "cmoveq 8(%[src]), %%rbx\n"
            "cmoveq 16(%[src]), %%r8\n"
            "cmovel %%edi, %%r9d\n"
            "movq %%rax, (%[dst])\n"
            "movq %%rbx, 8(%[dst])\n"
            "movq %%r8, 16(%[dst])\n"
            "movl %%r9d, (%[src])\n"
            "addq $24, %[src]\n"
            "addq $24, %[dst]\n"
            "movl $7, %%ecx\n"          // Loop 7 times moving 40 bytes at a time = 280 bytes
            ".stashExtractL%=:\n"
            "test %[flag], %[flag]\n"
            "movq (%[dst]), %%rax\n"
            "movq 8(%[dst]), %%rbx\n"
            "movq 16(%[dst]), %%rdi\n"
            "movq 24(%[dst]), %%r8\n"
            "movq 32(%[dst]), %%r9\n"
            "cmoveq (%[src]), %%rax\n"
            "cmoveq 8(%[src]), %%rbx\n"
            "cmoveq 16(%[src]), %%rdi\n"
            "cmoveq 24(%[src]), %%r8\n"
            "cmoveq 32(%[src]), %%r9\n"
            "movq %%rax, (%[dst])\n"
            "movq %%rbx, 8(%[dst])\n"
            "movq %%rdi, 16(%[dst])\n"
            "movq %%r8, 24(%[dst])\n"
            "movq %%r9, 32(%[dst])\n"
            "addq $40, %[dst]\n"
            "addq $40, %[src]\n"
            "dec %%ecx\n"
            "jne .stashExtractL%=\n"
            "inc %%ecx\n"               // ecx <- 1
            "test %[flag], %[flag]\n"
            "cmovel %%ecx, %[seen]\n"   // Move 1 into %[seen]
            : [seen] "+r" (seen), [dst] "+r" (dstblock), [src] "+r" (stash)
            : [flag] "r" ((stash->h.t == -1) | ((stash->h.leaf_id ^ leaf_id) & mask) | seen)
            : "cc", "memory", "rax", "rbx", "rcx", "rdi", "r8", "r9"
        );
    }
}

// An oblivious constant-time implementation of moving a block specified by
// a timestamp t.  If candidate->t == t, *candidate is copied to *dst and
// candidate->t is set to -1.
static inline void blockCondMove(Block *dst, Block *candidate, int t) {
    __asm__ volatile (
        "# If %[flag] is 1, move 304 bytes (38*8) from %[src] to %[dst]\n"
        "# and set *%[src] (== src->h.t) to -1\n"
        "xor %%edi, %%edi\n"
        "dec %%edi\n"               // edi <- -1
        "test %[flag], %[flag]\n"
        "movq (%[dst]), %%rax\n"    // Move the first 24 bytes and set *%[src] to -1
        "movq 8(%[dst]), %%rbx\n"
        "movq 16(%[dst]), %%r8\n"
        "cmovneq (%[src]), %%rax\n"
        "movl (%[src]), %%r9d\n"
        "cmovneq 8(%[src]), %%rbx\n"
        "cmovneq 16(%[src]), %%r8\n"
        "cmovnel %%edi, %%r9d\n"
        "movq %%rax, (%[dst])\n"
        "movq %%rbx, 8(%[dst])\n"
        "movq %%r8, 16(%[dst])\n"
        "movl %%r9d, (%[src])\n"
        "addq $24, %[src]\n"
        "addq $24, %[dst]\n"
        "movl $7, %%ecx\n"          // Loop 7 times moving 40 bytes at a time = 280 bytes
        ".blockCondMoveL%=:\n"
        "test %[flag], %[flag]\n"
        "movq (%[dst]), %%rax\n"
        "movq 8(%[dst]), %%rbx\n"
        "movq 16(%[dst]), %%rdi\n"
        "movq 24(%[dst]), %%r8\n"
        "movq 32(%[dst]), %%r9\n"
        "cmovneq (%[src]), %%rax\n"
        "cmovneq 8(%[src]), %%rbx\n"
        "cmovneq 16(%[src]), %%rdi\n"
        "cmovneq 24(%[src]), %%r8\n"
        "cmovneq 32(%[src]), %%r9\n"
        "movq %%rax, (%[dst])\n"
        "movq %%rbx, 8(%[dst])\n"
        "movq %%rdi, 16(%[dst])\n"
        "movq %%r8, 24(%[dst])\n"
        "movq %%r9, 32(%[dst])\n"
        "addq $40, %[dst]\n"
        "addq $40, %[src]\n"
        "dec %%ecx\n"
        "jne .blockCondMoveL%=\n"
        : [dst] "+r" (dst), [src] "+r" (candidate)
        : [flag] "r" (candidate->h.t == t)
        : "cc", "memory", "rax", "rbx", "rcx", "rdi", "r8", "r9"
    );
}

// If *candidate < *potential_min, copy *candidate into *potential_min (obliviously)
static inline void blockheader_min(BlockHeader *potential_min,
        const BlockHeader *candidate) {
    __asm__ (
        "# If [potmin]->t == -1, or ([cand]->t != -1 && [cand]->kt < [potmin]->kt)\n"
        "# then copy 16 bytes from *[cand] into *[potmin]\n"
        "# where kt is the 8-byte value ((k<<32) + t)\n"
        "movq (%[cand]), %%rax\n"
        "movq (%[potmin]), %%rbx\n"
        "movq 8(%[cand]), %%r8\n"
        "movq 8(%[potmin]), %%r9\n"
        "cmpq %%rbx, %%rax\n"
        "setl %%dl\n"             // dl <- [cand]->kt < [potmin]->kt
        "cmpl $-1, %%eax\n"
        "setne %%cl\n"            // cl <- ([cand]->t != -1)
        "and %%cl, %%dl\n"        // dl <- ([cand]->t != -1 && [cand]->kt < [potmin]->kt)
        "cmpl $-1, %%ebx\n"
        "sete %%cl\n"             // cl <- ([potmin]->t == -1)
        "or %%cl, %%dl\n"         // dl <- the final flag
        "test %%dl, %%dl\n"
        "cmovneq %%rax, %%rbx\n"
        "cmovneq %%r8, %%r9\n"
        "movq %%rbx, (%[potmin])\n"
        "movq %%r9, 8(%[potmin])\n"
        :
        : [potmin] "r" (potential_min), [cand] "r" (candidate)
        : "cc", "memory", "rax", "rbx", "rcx", "rdx", "r8", "r9"
    );
}

BlockHeader OHeapReadPathEviction::onlyFindMin() {
  BlockHeader potential_min = this->storage->readRoot().subtree_min;
  for (size_t i=0; i<this->short_stashsize;++i) {
    blockheader_min(&potential_min, &(this->stash[i].h));
  }
  return potential_min;
}

std::pair<int, int> OHeapReadPathEviction::onlyInsertBlock(int k, const unsigned char *v) {
  auto ret = std::pair<int, int>();

  int leaf = this->rand_gen->getRandomLeaf();
  if (stashInsert(this->stash.data(), this->short_stashsize, Block(leaf, this->t, k, v))) {
    throw std::runtime_error("Short stash full while inserting");
  }

  ret.first = leaf;
  ret.second = this->t;

  this->evictAfterInsert();

  return ret;
}

/** pedagogy */
void OHeapReadPathEviction::insertDummyBlock() {
  this->evictAfterInsert();
}

/** convenience function */
void OHeapReadPathEviction::evictAfterInsert() {
  // two random paths non-overlapping except at the root
  // = random leaf_id in the first half and random leaf_id in the second half (over all leaf_ids)
  int random_leaf_one = this->rand_gen->getRandomLeaf()/2;
  int random_leaf_two = this->rand_gen->getRandomLeaf()/2 + this->num_leaves/2;

  this->evictAndUpdateMin(random_leaf_one);
  this->evictAndUpdateMin(random_leaf_two);
}

Block OHeapReadPathEviction::onlyDeleteBlock(int leaf, long t) {
  // try to delete from storage
  Block deleted_block = this->readAndRemove(leaf, t);

  // try to delete from stash
  Block *stashptr = this->stash.data();
  for (size_t i = 0; i < this->short_stashsize; i++) {
    blockCondMove(&deleted_block, stashptr, t);
    ++stashptr;
  }

  return deleted_block;
}

void OHeapReadPathEviction::evictAndUpdateMin(int leaf) {
  if (leaf == -1) {
    // This can only happen when type_hiding_security is set
    leaf = this->rand_gen->getRandomLeaf();
  }
  auto paths = this->storage->readPath(leaf);
  this->evictAndUpdateMin(leaf, paths.first, paths.second);
  this->storage->writebackPath();
}

void OHeapReadPathEviction::evictAndUpdateMin(int leaf, std::vector<Bucket>& path,
        const std::vector<Bucket>& pathsiblings) {
  // pull all blocks on path into stash
  // At this point, we know that at most the first short_stashsize elements
  // of the stash are occupied.  So until we hit the long_stashsize, we
  // can just memmove the blocks in the bucket into the stash.  Once we
  // hit the long_stashsize, we'll need to obliviously insert the rest of
  // the blocks.
  Block *stashstart = this->stash.data();
  size_t stashoffset = this->short_stashsize;
  for (int l = 0; l < num_levels; l++) {
    const Block* blocks = path[l].blocks;
    if (stashoffset + Bucket::BUCKET_CAPACITY <= this->long_stashsize) {
        memmove(stashstart+stashoffset, blocks, Bucket::BUCKET_CAPACITY * sizeof(Block));
    } else {
        for (int i = 0; i < Bucket::BUCKET_CAPACITY; ++i) {
            if (stashInsert(stashstart, this->long_stashsize, blocks[i])) {
              throw std::runtime_error("Long stash full while evicting");
            }
        }
    }
    stashoffset += Bucket::BUCKET_CAPACITY;
  }

  // push as many stash blocks as possible into the leaf bucket
  // then push as many stash blocks as possible into the leaf's parent
  // etc.
  unsigned int mask = 0xffffffff;
  for (int l = num_levels-1; l >= 0; l--) {
    Bucket *bucket = &path[l];
    for (int i = 0; i < Bucket::BUCKET_CAPACITY; ++i) {
        stashExtract(stashstart, this->long_stashsize, leaf, mask,
            bucket->blocks + i);
    }
    // Remove the lowest set bit from the mask
    mask &= (mask-1);
  }

  // Obliviously compress the long stash into the short stash
  Block *stashend = stashstart+(this->long_stashsize);
  for (Block *block = stashstart+(this->short_stashsize); block < stashend; ++block) {
    if (stashInsert(stashstart, this->short_stashsize, *block)) {
      throw std::runtime_error("Short stash full while compressing");
    }
    block->h.t = -1;
  }

  this->updateMin(leaf, path, pathsiblings);
}

Block OHeapReadPathEviction::readAndRemove(int leaf, long t) {
    int randleaf = rand_gen->getRandomLeaf();
    // if (leaf == -1) { leaf = randleaf; }, but obliviously
    __asm__(
        "test %[flag], %[flag]\n"
        "cmovnel %[randleaf], %[leaf]\n"
        :
        : [randleaf] "r" (randleaf), [leaf] "r" (leaf), [flag] "r" (leaf == -1)
        : "cc"
    );

    Block ret;
    auto paths = this->storage->readPath(leaf);
    for (int l = 0; l < this->num_levels; l++) {
        Block *blocks = paths.first[l].blocks;
        for (int i = 0; i < Bucket::BUCKET_CAPACITY; ++i) {
            blockCondMove(&ret, blocks+i, t);
        }
    }
    this->evictAndUpdateMin(leaf, paths.first, paths.second);
    this->storage->writebackPath();
    return ret;
}

void OHeapReadPathEviction::updateMin(int leaf, std::vector<Bucket>& path,
        const std::vector<Bucket>& pathsiblings) {
    for (int l = num_levels-1; l >= 0; l--) {
      BlockHeader *potential_min = &(path[l].subtree_min);
      potential_min->t = -1;

      const Block* blocks_on_path = path[l].blocks;
      for (int i=0;i<Bucket::BUCKET_CAPACITY; ++i) {
        blockheader_min(potential_min, &(blocks_on_path[i].h));
      }
      if (l != num_levels-1) {
        blockheader_min(potential_min, &path[l+1].subtree_min);
        blockheader_min(potential_min, &pathsiblings[l].subtree_min);
      }
    }
}

int OHeapReadPathEviction::P(int leaf, int level) {
  if (level == 0) return 0;

  return (1<<level)-1 + (leaf>>(num_levels-1-level));
}

std::vector<Block>& OHeapReadPathEviction::getStash() {
  return this->stash;
}

int OHeapReadPathEviction::getStashSize() {
  return this->stash.size();
}

int OHeapReadPathEviction::getNumLeaves() {
  return this->num_leaves;
}

int OHeapReadPathEviction::getNumLevels() {
  return this->num_levels;
}

int OHeapReadPathEviction::getNumBlocks() {
  return this->num_blocks;
}

int OHeapReadPathEviction::getNumBuckets() {
  return this->num_buckets;
}

/*
std::ostream& operator<<(std::ostream& out, const OHeapReadPathEviction& oheap) {
  out << "==========" << std::endl;
  for (int i = 0; i < oheap.num_buckets; i++) {
    std::vector<Block> blocks = oheap.storage->readBucket(i).getBlocks();
    bool all_dummy = true;
    for (const Block& bl : blocks) {
      if (bl.leaf_id != -1) {
        all_dummy = false;
        break;
      }
    }
    if (!all_dummy) {
      out << "Bucket " << i << " contains:" << std::endl;

      for (int j = 0; j < oheap.bucket_size; j++) {
        const Block& bl = blocks.at(j);
        if (bl.leaf_id != -1) {
          out << bl << std::endl;
        }
      }

      out << std::endl;
    }
  }
  out << "==========" << std::endl;
  return out;
}
*/
