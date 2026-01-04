#ifndef POHEAP_EXTERNALSERVERSTORAGE_H
#define POHEAP_EXTERNALSERVERSTORAGE_H

#include "UntrustedStorageInterface.h"
#include <sgx_tcrypto.h>

class ExternalServerStorage : public UntrustedStorageInterface {
public:
    ExternalServerStorage();
    ~ExternalServerStorage();

    /** Must be called before reading or writing any buckets. Cannot be changed after being set. */
    void setCapacity(int total_num_of_buckets);

    /** Returns a reference to the root of the tree. */
    Bucket& readRoot();

    /** Reads the path to the given leaf into the local path, and returns a reference to that path. */
    std::pair<std::vector<Bucket>&,const std::vector<Bucket>&> readPath(int leaf);

    /** Write back the local path into storage, invalidating the local path. */
    void writebackPath();

private:
    int capacity;
    int num_levels;
    int num_leaves;
    sgx_aes_gcm_128bit_key_t key;
    unsigned char *untrusted_encbuckets;
    unsigned char *untrusted_ivs;
    sgx_aes_gcm_128bit_tag_t *untrusted_tags;
    Bucket root;  // The root bucket is just kept in enclave memory, since it's accessed a lot
    sgx_aes_gcm_128bit_tag_t lefttag, righttag;  // The GCM tags of the children of the root
    std::vector<Bucket> curpath;
    std::vector<Bucket> pathsiblings;
    std::vector<unsigned char> path_and_sibling_tags;
    int curleaf;

    void init_subtree(int index, int level, const Bucket &empty_bucket,
            sgx_aes_gcm_128bit_tag_t *subtree_root_tag);
};

#endif
