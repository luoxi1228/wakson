#include "ServerStorage.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include "utils.hpp"

bool ServerStorage::is_initialized = false;
bool ServerStorage::is_capacity_set = false;

ServerStorage::ServerStorage() : buckets() {
  if (ServerStorage::is_initialized) {
    throw std::runtime_error("Only one ServerStorage can be used at a time in this implementation!");
  }
  ServerStorage::is_initialized = true;
}

void ServerStorage::setCapacity(int total_num_of_buckets) {
  if (this->is_capacity_set) {
    throw std::runtime_error("Capacity of ServerStorage cannot be changed");
  }
  this->is_capacity_set = true;
  this->capacity = total_num_of_buckets;
  this->buckets.assign(total_num_of_buckets, Bucket());
  this->num_levels = ceil(log2(total_num_of_buckets));
  this->num_leaves = 1<<(this->num_levels - 1);
  this->curleaf = -1;
  this->curpath.assign(this->num_levels, Bucket());
  this->pathsiblings.assign(this->num_levels-1, Bucket());
}

void ServerStorage::checkRWValid(int position) {
  if (!this->is_capacity_set) {
    throw std::runtime_error("Please call setCapacity before reading or writing any block");
  }
  if (position >= this->capacity || position < 0) {
    throw std::runtime_error("You are trying to access Bucket " + std::to_string(position) +
        ", but this Server contains only " + std::to_string(this->capacity) + " buckets.");
  }
}

Bucket& ServerStorage::readRoot() {
  return this->buckets[0];
}

Bucket ServerStorage::readBucket(int position) {
  this->checkRWValid(position);
  return this->buckets.at(position);
}

void ServerStorage::writeBucket(int position, const Bucket& bucket_to_write) {
  this->checkRWValid(position);
  this->buckets.at(position) = bucket_to_write;
}

std::pair<std::vector<Bucket>&,const std::vector<Bucket>&>
        ServerStorage::readPath(int leaf) {
  if (this->curleaf >= 0) {
    printf("Trying to read multiple paths at once\n");
  }
  // idx is 1-based (the root is 1, its children are 2 and 3, 2's children are 4 and 5, etc.
  // the benefit of this notation is that the parent of idx is just int(idx/2), and the sibling
  // of idx is idx^1.
  int idx = (1<<(this->num_levels-1))+leaf;
  for (int level = this->num_levels-1; level>=0; --level) {
    this->curpath[level] = this->buckets[idx-1];
    if (level > 0) {
        this->pathsiblings[level-1] = this->buckets[(idx^1)-1];
    }
    idx >>= 1;
  }
  this->curleaf = leaf;
  return std::pair<std::vector<Bucket>&,const std::vector<Bucket>&>(
    this->curpath, this->pathsiblings);
}

void ServerStorage::writebackPath() {
  if (this->curleaf < 0) {
    printf("No current path to write back");
  }
  // idx is 1-based (the root is 1, its children are 2 and 3, 2's children are 4 and 5, etc.
  // the benefit of this notation is that the parent of idx is just int(idx/2)
  int idx = (1<<(this->num_levels-1))+this->curleaf;
  for (int level = this->num_levels-1; level>=0; --level) {
    this->buckets[idx-1] = this->curpath[level];
    idx >>= 1;
  }
  this->curleaf = -1;
}
