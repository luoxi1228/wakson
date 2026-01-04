#include <array>
#include <sgx_tcrypto.h>
#include "HeapSort.h"
#include "ServerStorage.h"
#include "ExternalServerStorage.h"
#include "SGXRand.h"
#include "OHeapReadPathEviction.h"
#include "Block.h"
#include "utils.hpp"

static unsigned char aesdeckey[SGX_AESGCM_KEY_SIZE];
static unsigned char aesenckey[SGX_AESGCM_KEY_SIZE];

// Load test AES keys into the enclave
void HeapSort_LoadTestKeys(unsigned char inkey[16], unsigned char outkey[16]) {
    memmove(aesdeckey, inkey, 16);
    memmove(aesenckey, outkey, 16);
}

// The output buffer and input buffer are each of length len bytes, which must
// be a multiple of 320.  The input buffer should contain 0 or more 320-byte
// chunks, each of which is a GCM encryption: 12 bytes IV, 288 bytes data, 16
// bytes tag.  The 288-byte data chunks, when decrypted, consist of a 4-byte
// key followed by a 288-byte data block.  The output buffer will be filled
// with reencryptions the same chunks, sorted by key.  The internal flag
// determines whether block storage will be internal or external to enclave
// memory.
void HeapSort(unsigned char *outbuf, const unsigned char *inbuf, size_t len, int internal) {
    long t0, t1, t2, t3;
    ocall_clock(&t0);
    size_t chunk_size = (sizeof(int) + Block::BLOCK_SIZE);
    size_t enc_chunk_size = SGX_AESGCM_IV_SIZE + chunk_size + SGX_AESGCM_MAC_SIZE;
    int num_blocks = len / enc_chunk_size;
    UntrustedStorageInterface* storage;
    if (internal) {
        storage = new ServerStorage();
    } else {
        storage = new ExternalServerStorage();
    }
    RandForOHeapInterface* random = new SGXRand();
    OHeapInterface *oheap = new OHeapReadPathEviction(storage, random, num_blocks, false);

    try {
        ocall_clock(&t1);

        // Insert the blocks into the oheap
        const unsigned char *inbufend = inbuf + num_blocks * enc_chunk_size;
        for (const unsigned char *encchunkptr = inbuf; encchunkptr < inbufend; encchunkptr += enc_chunk_size) {
            unsigned char chunkdata[chunk_size];
            sgx_status_t aesret = sgx_rijndael128GCM_decrypt(
                &aesdeckey, encchunkptr + SGX_AESGCM_IV_SIZE, chunk_size,
                chunkdata, encchunkptr, SGX_AESGCM_IV_SIZE, NULL, 0,
                (const sgx_aes_gcm_128bit_tag_t*)(encchunkptr + SGX_AESGCM_IV_SIZE + chunk_size));
            if (aesret != SGX_SUCCESS) {
                printf("sgx_rijndael128GCM_decrypt failure (%x)\n", aesret);
                return;
            }
            int key = *(int *)chunkdata;
            oheap->insertBlock(key, chunkdata+sizeof(int));
        }

        ocall_clock(&t2);

        // Extract the sorted blocks from the oheap
        unsigned char *outbufend = outbuf + num_blocks * enc_chunk_size;
        for (unsigned char *encchunkptr = outbuf; encchunkptr < outbufend; encchunkptr += enc_chunk_size) {
            Block b = oheap->extractMin();
            unsigned char chunkdata[chunk_size];
            *(int *)chunkdata = b.h.k;
            memmove(chunkdata+sizeof(int), b.v, Block::BLOCK_SIZE);
            sgx_status_t randret = sgx_read_rand(encchunkptr, SGX_AESGCM_IV_SIZE);
            if (randret != SGX_SUCCESS) {
                printf("sgx_read_rand failure (%x)\n", randret);
                return;
            }
            sgx_status_t aesret = sgx_rijndael128GCM_encrypt(
                &aesenckey, chunkdata, chunk_size,
                encchunkptr + SGX_AESGCM_IV_SIZE,
                encchunkptr, SGX_AESGCM_IV_SIZE, NULL, 0,
                (sgx_aes_gcm_128bit_tag_t*)(encchunkptr + SGX_AESGCM_IV_SIZE + chunk_size));
            if (aesret != SGX_SUCCESS) {
                printf("sgx_rijndael128GCM_encrypt failure (%x)\n", aesret);
                return;
            }
        }

        ocall_clock(&t3);

        // CLOCKS_PER_SEC == 1000000, so CLOCKS_PER_MS == 1000
        double setup_ms = ((double)(t1-t0))/1000.0;
        double insert_ms = ((double)(t2-t1))/1000.0;
        double extract_ms = ((double)(t3-t2))/1000.0;

        printf("Setup time %lf ms\nInsert time %lf ms\nExtract time %lf ms\n", setup_ms, insert_ms, extract_ms);
    } catch (std::runtime_error e) {
        printf("Exception: %s\n", e.what());
    }

    delete oheap;
    delete random;
    delete storage;

}
