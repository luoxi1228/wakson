#include "ExternalServerStorage.h"
#include <cmath>
#include <stdexcept>
#include "utils.hpp"

ExternalServerStorage::ExternalServerStorage() :
        capacity(0), num_levels(0), num_leaves(0),
        untrusted_encbuckets(NULL), untrusted_ivs(NULL),
        untrusted_tags(NULL), root(), curpath(), curleaf(-1) {
}

void ExternalServerStorage::setCapacity(int total_num_of_buckets) {
    this->capacity = total_num_of_buckets;
    this->num_levels = ceil(log2(total_num_of_buckets));
    this->num_leaves = 1<<(this->num_levels - 1);
    this->curleaf = -1;
    this->curpath.assign(this->num_levels, Bucket());
    this->pathsiblings.assign(this->num_levels-1, Bucket());
    this->path_and_sibling_tags.reserve(2*(this->num_levels-1)*SGX_AESGCM_MAC_SIZE);
    int num_external_buckets = (1<<(this->num_levels)) - 2;  // The root is stored internally

    // Allocate external storage
    untrustedMemAllocate(&this->untrusted_encbuckets,
            num_external_buckets * sizeof(Bucket));
    untrustedMemAllocate(&this->untrusted_ivs,
            num_external_buckets * SGX_AESGCM_IV_SIZE);
    untrustedMemAllocate((unsigned char **)(&this->untrusted_tags),
            num_external_buckets * sizeof(sgx_aes_gcm_128bit_tag_t));

    // Pick a random AES GCM key for external storage
    sgx_read_rand((unsigned char *)(&this->key), sizeof(this->key));

    // Fill the tree with empty buckets
    Bucket empty_bucket;
    init_subtree(2, 1, empty_bucket, &(this->lefttag));
    init_subtree(3, 1, empty_bucket, &(this->righttag));
}

// Indices are 1-based: the root is 1, its children are 2 and 3, 2's children
// are 4 and 5, etc.  The parent of idx is floor(idx/2); the children of idx are
// (2*idx) and (2*idx+1); the sibling of idx is (idx^1).
// Levels are 0-based: the root is at level 0, its children are at level 1
void ExternalServerStorage::init_subtree(int index, int level,
        const Bucket &empty_bucket,
        sgx_aes_gcm_128bit_tag_t *subtree_root_tag) {
    sgx_aes_gcm_128bit_tag_t childtags[2];
    unsigned char iv[SGX_AESGCM_IV_SIZE];
    // Post-order traversal; we include our children's GCM tags as our
    // own GCM associated data
    sgx_read_rand(iv, sizeof(iv));
    memmove(this->untrusted_ivs + (index-2)*SGX_AESGCM_IV_SIZE, iv,
            SGX_AESGCM_IV_SIZE);
    if (level < this->num_levels - 1) {
        init_subtree(2*index, level+1, empty_bucket, &(childtags[0]));
        init_subtree(2*index+1, level+1, empty_bucket, &(childtags[1]));
        sgx_rijndael128GCM_encrypt(&this->key,
                (const unsigned char*)&empty_bucket, sizeof(Bucket),
                this->untrusted_encbuckets + (index-2)*sizeof(Bucket),
                iv, SGX_AESGCM_IV_SIZE, (unsigned char *)childtags,
                sizeof(childtags), subtree_root_tag);
    } else {
        // No associated data for leaf nodes
        sgx_rijndael128GCM_encrypt(&this->key,
                (const unsigned char*)&empty_bucket, sizeof(Bucket),
                this->untrusted_encbuckets + (index-2)*sizeof(Bucket),
                iv, SGX_AESGCM_IV_SIZE, NULL, 0, subtree_root_tag);
    }
    memmove(this->untrusted_tags + (index-2),
        (unsigned char *)subtree_root_tag, sizeof(sgx_aes_gcm_128bit_tag_t));
}

ExternalServerStorage::~ExternalServerStorage() {
    untrustedMemFree(this->untrusted_encbuckets);
    untrustedMemFree(this->untrusted_ivs);
    untrustedMemFree((unsigned char *)(this->untrusted_tags));
}

Bucket& ExternalServerStorage::readRoot() {
    return this->root;
}

std::pair<std::vector<Bucket>&,const std::vector<Bucket>&>
        ExternalServerStorage::readPath(int leaf) {
    // Read the path (except the root) and its siblings into path and
    // pathsiblings
    unsigned char *localtagdata = path_and_sibling_tags.data();
    Bucket *curpathdata = curpath.data();
    Bucket *pathsiblingsdata = pathsiblings.data();

    // idx is 1-based (the root is 1, its children are 2 and 3, 2's children
    // are 4 and 5, etc.  the benefit of this notation is that the parent of
    // idx is just int(idx/2), and the sibling of idx is idx^1.
    int idx = (1<<(this->num_levels-1))+leaf;
    for (int l = this->num_levels-1; l >= 1; --l) {
        // Read bucket idx into the path, and its sibling into pathsiblings
        int siblingidx = idx^1;
        // Read the purported GCM tags into enclave memory
        sgx_aes_gcm_128bit_tag_t *nodetag =
                (sgx_aes_gcm_128bit_tag_t*)(localtagdata) +
                    (2*(l-1) + (idx&1));
        sgx_aes_gcm_128bit_tag_t *siblingtag =
                (sgx_aes_gcm_128bit_tag_t*)(localtagdata) +
                    (2*(l-1) + (siblingidx&1));
        memmove(nodetag, untrusted_tags + (idx-2), SGX_AESGCM_MAC_SIZE);
        memmove(siblingtag, untrusted_tags + (siblingidx-2), SGX_AESGCM_MAC_SIZE);

        const unsigned char *nodeiv =
                untrusted_ivs + (idx-2)*SGX_AESGCM_IV_SIZE;
        const unsigned char *siblingiv =
                untrusted_ivs + (siblingidx-2)*SGX_AESGCM_IV_SIZE;

        const unsigned char *nodeenc =
                untrusted_encbuckets + (idx-2)*sizeof(Bucket);
        const unsigned char *siblingenc =
                untrusted_encbuckets + (siblingidx-2)*sizeof(Bucket);

        sgx_status_t aesret, sibaesret;
        if (l < this->num_levels - 1) {
            // The AAD is the concatenation of the left and right child's
            // GCM tags
            const unsigned char *aad = localtagdata +
                2*l*SGX_AESGCM_MAC_SIZE;
            const unsigned char *siblingaad =
                (const unsigned char *)(untrusted_tags + (2*siblingidx-2));
            aesret = sgx_rijndael128GCM_decrypt(
                &this->key, nodeenc, sizeof(Bucket),
                (unsigned char *)(curpathdata+l),
                nodeiv, SGX_AESGCM_IV_SIZE,
                aad, SGX_AESGCM_MAC_SIZE*2, nodetag);
            sibaesret = sgx_rijndael128GCM_decrypt(
                &this->key, siblingenc, sizeof(Bucket),
                (unsigned char *)(pathsiblingsdata+(l-1)),
                siblingiv, SGX_AESGCM_IV_SIZE,
                siblingaad, SGX_AESGCM_MAC_SIZE*2, siblingtag);
        } else {
            // This is a leaf node; the AAD is empty
            aesret = sgx_rijndael128GCM_decrypt(
                &this->key, nodeenc, sizeof(Bucket),
                (unsigned char *)(this->curpath.data()+l),
                nodeiv, SGX_AESGCM_IV_SIZE,
                NULL, 0, nodetag);
            sibaesret = sgx_rijndael128GCM_decrypt(
                &this->key, siblingenc, sizeof(Bucket),
                (unsigned char *)(this->pathsiblings.data()+(l-1)),
                siblingiv, SGX_AESGCM_IV_SIZE,
                NULL, 0, siblingtag);
        }
        if (aesret != SGX_SUCCESS) {
            printf("failed at %d %d\n", l, idx);
            for(int i=0;i<12;++i) printf("%02x", nodeiv[i]); printf(" <- %p\n", nodeiv);
            for(int i=0;i<sizeof(Bucket);++i) printf("%02x", nodeenc[i]); printf("\n");
            for(int i=0;i<SGX_AESGCM_MAC_SIZE*2;++i) printf("%02x", (localtagdata + 2*l*SGX_AESGCM_MAC_SIZE)[i]); printf(" <- %p\n", nodetag);
            throw std::runtime_error("Node decrypt failure");
        }
        if (sibaesret != SGX_SUCCESS) {
            printf("failed at %d %d\n", l, idx);
            for(int i=0;i<12;++i) printf("%02x", siblingiv[i]); printf(" <- %p\n", siblingiv);
            for(int i=0;i<sizeof(Bucket);++i) printf("%02x", siblingenc[i]); printf("\n");
            for(int i=0;i<SGX_AESGCM_MAC_SIZE*2;++i) printf("%02x", ((unsigned char *)(untrusted_tags + (2*siblingidx-2)))[i]); printf(" <- %p\n", siblingtag);
            throw std::runtime_error("Sibling decrypt failure");
        }
        idx >>= 1;
    }
    // Copy the root into the path
    this->curpath[0] = this->root;
    // Check the top-level GCM tags
    if (memcmp(this->lefttag, localtagdata, SGX_AESGCM_MAC_SIZE) ||
        memcmp(this->righttag, localtagdata+SGX_AESGCM_MAC_SIZE,
                SGX_AESGCM_MAC_SIZE)) {
        throw std::runtime_error("Integrity check failed");
    }
    this->curleaf = leaf;
    return std::pair<std::vector<Bucket>&,const std::vector<Bucket>&>(
        this->curpath, this->pathsiblings);
}

void ExternalServerStorage::writebackPath() {
    unsigned char *localtagdata = path_and_sibling_tags.data();
    Bucket *curpathdata = curpath.data();
    int idx = (1<<(this->num_levels-1))+this->curleaf;
    for (int l = this->num_levels-1; l >= 1; --l) {
        int siblingidx = idx^1;
        unsigned char iv[SGX_AESGCM_IV_SIZE];

        sgx_aes_gcm_128bit_tag_t *nodetag =
                (sgx_aes_gcm_128bit_tag_t*)(localtagdata) +
                    (2*(l-1) + (idx&1));
        unsigned char *extnodetag =
            (unsigned char *)(untrusted_tags + (idx-2));

        const unsigned char *nodeiv =
                untrusted_ivs + (idx-2)*SGX_AESGCM_IV_SIZE;

        const unsigned char *nodeenc =
                untrusted_encbuckets + (idx-2)*sizeof(Bucket);

        sgx_read_rand(iv, sizeof(iv));
        memmove(this->untrusted_ivs + (idx-2)*SGX_AESGCM_IV_SIZE, iv,
                SGX_AESGCM_IV_SIZE);

        if (l < this->num_levels - 1) {
            // A non-leaf node; the AAD is the concatenation of the children's
            // GCM tags
            const unsigned char *aad = localtagdata +
                2*l*SGX_AESGCM_MAC_SIZE;
            sgx_rijndael128GCM_encrypt(&this->key,
                (unsigned char*)(curpathdata+l), sizeof(Bucket),
                this->untrusted_encbuckets + (idx-2)*sizeof(Bucket),
                iv, SGX_AESGCM_IV_SIZE, aad, 2*SGX_AESGCM_MAC_SIZE,
                nodetag);
        } else {
            // This is a leaf node; no AAD
            sgx_rijndael128GCM_encrypt(&this->key,
                (unsigned char*)(curpathdata+l), sizeof(Bucket),
                this->untrusted_encbuckets + (idx-2)*sizeof(Bucket),
                iv, SGX_AESGCM_IV_SIZE, NULL, 0, nodetag);
        }
        memmove(extnodetag, nodetag, SGX_AESGCM_MAC_SIZE);

        idx >>= 1;
    }
    // Store the new root and its childrens' GCM tags locally
    this->root = this->curpath[0];
    memmove(this->lefttag, localtagdata, SGX_AESGCM_MAC_SIZE);
    memmove(this->righttag, localtagdata+SGX_AESGCM_MAC_SIZE,
            SGX_AESGCM_MAC_SIZE);
    this->curleaf = -1;
}
