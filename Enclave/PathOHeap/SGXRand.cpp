#include <sgx_trts.h>
#include <stdexcept>
#include "SGXRand.h"
#include "utils.hpp"

SGXRand::SGXRand() {
}

void SGXRand::setBound(int num_leaves) {
    unsigned int maxval;
    if (num_leaves < 1) {
	num_leaves = 1;
    }
    maxval = num_leaves - 1;

    // Find the smallest number of the form 2^n - 1 that is at least maxval
    unsigned int mask = 0;
    while (mask < maxval) {
        mask = (2*mask) + 1;
    }
    this->maxval = maxval;
    this->mask = mask;
    // printf("maxval = %d, mask = %08x\n", this->maxval, mask);
}

int SGXRand::getRandomLeaf() {
    unsigned char randbuf[4];
    while(1) {
        sgx_status_t res = sgx_read_rand(randbuf, 4);
        if (res != SGX_SUCCESS) {
            throw std::runtime_error("sgx_read_rand failed");
        }
        unsigned int randval = *(unsigned int*)randbuf;
        randval &= this->mask;
        if (randval <= this->maxval) {
            return int(randval);
        }
    }
}

void TestRand() {
    SGXRand r;

    for (int numleaves = -2; numleaves < 70; ++numleaves) {
        r.setBound(numleaves);
        int numrands = numleaves;
        if (numrands < 5) numrands = 5;
        for (int i = 0; i < numrands; ++i) {
            int rval = r.getRandomLeaf();
            printf("rand = %d\n", rval);
        }
    }
}
