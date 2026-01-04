#ifndef __HEAPSORT_H__
#define __HEAPSORT_H__

#ifdef __cplusplus
extern "C" {
#endif

// The output buffer and input buffer are each of length len bytes, which must
// be a multiple of 320.  The input buffer should contain 0 or more 320-byte
// chunks, each of which is a GCM encryption: 12 bytes IV, 288 bytes data, 16
// bytes tag.  The 288-byte data chunks, when decrypted, consist of a 4-byte
// key followed by a 288-byte data block.  The output buffer will be filled
// with reencryptions the same chunks, sorted by key.  The internal flag
// determines whether block storage will be internal or external to enclave
// memory.
void HeapSort(unsigned char *outbuf, const unsigned char *inbuf, size_t len, int internal);

// Load test AES keys into the enclave
void HeapSort_LoadTestKeys(unsigned char inkey[16], unsigned char outkey[16]);

#ifdef __cplusplus
}
#endif

#endif
