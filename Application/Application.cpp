#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <openssl/err.h>
#include "Application.hpp"
#include "gcm.h"

int main(int argc, char *argv[]){

  // Initialize libcrypto
  OpenSSL_add_all_algorithms();
  ERR_load_crypto_strings();

  clock_t process_start, process_stop;
  double total_time=0, ptime;
  OLib_initialize();
 
  /*
  int ms_node_id =  OSort_initializeMergeSplitNode(PARAM_F, PARAM_B, 0);

  for(int i=0; i<PARAM_N; i++){ 
    packet *packet_to_process = generate_packet(PARAM_L, PARAM_K, PARAM_F);
    unsigned char *packet_in_out = serialize_packet(packet_to_process);
    int returned_stream;

    #ifdef PRINT_PACKETS
      printf("Sending Packet:\n");
      displayPacket(packet_to_process);  
    #endif

    process_start = clock();
    returned_stream = OSort_processPacket(ms_node_id, packet_in_out, SERIALIZED_PACKET_SIZE);
    process_stop = clock();
    ptime = double(process_stop-process_start)/double(CLOCKS_PER_MS);
    total_time+=ptime;

    packet *packet_recieved = deserialize_packet(packet_in_out);

    #ifdef PRINT_PACKETS
      printf("Out Packet after processing:\n");
      displayPacket(packet_recieved);
      printf("\n\n");
    #endif

    if(packet_recieved->out_address!=returned_stream&&packet_recieved->out_address!=UINT64_MAX)
      printf("Failed: Received packet's output address and returned stream mismatch.\n");

    free(packet_to_process);
    free(packet_recieved); 
    free(packet_in_out);
  }

  printf("Total Time = %f ms\n", total_time);
  */

  // Test the SGX Enclave RNG
  //PathOHeap_TestRand();

  // Test HeapSort
  const size_t BLOCK_SIZE = 288;
  const size_t num_chunks = 10000;

  int randfd = open("/dev/urandom", O_RDONLY);
  if (randfd < 0) {
    throw std::runtime_error("Cannot open /dev/urandom");
  }

  // Load test AES keys into the enclave
  unsigned char inkey[16];
  unsigned char outkey[16];
  read(randfd, inkey, 16);
  read(randfd, outkey, 16);
  PathOHeap_HeapSort_LoadTestKeys(inkey, outkey);

  // An AES key for constructing and verifying the _plaintext_ of the blocks
  // so that we can ensure that the blocks come out unaltered
  unsigned char datakey[16];
  read(randfd, datakey, 16);
  // The plaintext data will itself be an AES-GCM encryption of all zeros: 12
  // bytes IV, BLOCK_SIZE-12-16 bytes ciphertext, 16 bytes tag
  unsigned char zeroes[BLOCK_SIZE-12-16];
  memset(zeroes, 0, sizeof(zeroes));

  const size_t chunk_size = sizeof(int) + BLOCK_SIZE;
  const size_t enc_chunk_size = 12 + chunk_size + 16;
  unsigned char inoutbuf[num_chunks * enc_chunk_size];
  unsigned char *bufend = inoutbuf + num_chunks * enc_chunk_size;
  for (unsigned char *encchunkptr = inoutbuf; encchunkptr < bufend; encchunkptr += enc_chunk_size) {
    unsigned char chunkbuf[chunk_size];
    // Pick a random chunk key
    read(randfd, chunkbuf, 4);
    // Set the plaintext of the message as above so that it's something we can
    // check was unaltered when it comes back out.  The plaintext will itself be
    // an AES-GCM ciphertext
    // Choose a random IV
    read(randfd, chunkbuf+4, 12);
    // GCM-encrypt, using the chunk key as the associated data
    if (sizeof(zeroes) != gcm_encrypt(zeroes, sizeof(zeroes), chunkbuf, 4, datakey,
            chunkbuf+4, 12, chunkbuf+4+12, chunkbuf+4+12+sizeof(zeroes))) {
        printf("Inner encryption failed\n");
    }
    // Encrypt the chunk to the enclave
    read(randfd, encchunkptr, 12);
    if (chunk_size != gcm_encrypt(chunkbuf, chunk_size, NULL, 0, inkey,
            encchunkptr, 12, encchunkptr+12, encchunkptr+12+chunk_size)) {
        printf("Encryption failed\n");
        break;
    }
  }
  close(randfd);

  // Sort it
  const char *use_internal_storage_env = getenv("PATHOHEAP_INTSTORAGE");
  int use_internal_storage = 0;
  if (use_internal_storage_env) {
    use_internal_storage = atoi(use_internal_storage_env);
  }
  process_start = clock();
  PathOHeap_HeapSort(inoutbuf, inoutbuf, num_chunks * enc_chunk_size, use_internal_storage);
  process_stop = clock();
  ptime = double(process_stop-process_start)/double(CLOCKS_PER_MS);

  #if 1
  // Decrypt the output and do a sanity check
  int lastkey = INT_MIN;
  int numfailed = 0;
  int cnum = 0;
  for (unsigned char *encchunkptr = inoutbuf; encchunkptr < bufend; encchunkptr += enc_chunk_size) {
    unsigned char chunkbuf[chunk_size];
    ++cnum;
    if (chunk_size != gcm_decrypt(encchunkptr+12, chunk_size, NULL, 0, encchunkptr+12+chunk_size, outkey, encchunkptr, 12, chunkbuf)) {
        printf("Decryption failed %d/%d\n", ++numfailed, cnum);
    }

    int chunkkey = *(int *)chunkbuf;
    if (chunkkey < lastkey) {
        printf("Unsorted!\n");
    }
    // Check the integrity of the data that came out
    unsigned char should_be_zeroes[sizeof(zeroes)];
    if (sizeof(zeroes) != gcm_decrypt(chunkbuf+4+12, sizeof(zeroes), chunkbuf, 4,
            chunkbuf+4+12+sizeof(zeroes), datakey, chunkbuf+4, 12, should_be_zeroes) ||
            memcmp(should_be_zeroes, zeroes, sizeof(zeroes))) {
        printf("Inner decryption failed\n");
        for(int i=0;i<4;++i) printf("%02x", chunkbuf[i]); printf("\n");
        for(int i=0;i<12;++i) printf("%02x", chunkbuf[4+i]); printf("\n");
        for(int i=0;i<sizeof(zeroes);++i) printf("%02x", chunkbuf[4+12+i]); printf("\n");
        for(int i=0;i<16;++i) printf("%02x", chunkbuf[4+12+sizeof(zeroes)+i]); printf("\n");
    }
    lastkey = chunkkey;
  }
  #endif
  printf("Sort time = %f ms\n", ptime);

  /*
  //Test Bucket Oblivious Sort

  //TODO: For now set params into this params structure and pass it to 
  //      BucketOSort_setParams()
  bos_params params;
  params.n = PARAM_N;
  params.f = PARAM_F;
  // NOTE: b = 2n/z is the number of buckets 
  // but we also want b to be a perfect power of f
  // For now we select n such that this is true!
  params.b = ((2 * PARAM_N)/PARAM_Z); 
  params.Z = PARAM_Z;
  // b = f^(d-1) or ie d = log_f(b) - 1
  params.d = PARAM_D; 
  params.B = PARAM_B;

  //TODO: Generate PARAM_N packets, serialized and push them into inbuf
  uint64_t buf_len = (uint64_t)(PARAM_N) * SERIALIZED_PACKET_SIZE; 
  printf("buf_len = %ld\n", buf_len);
  unsigned char *inbuf = (unsigned char*) malloc(buf_len);
  if(inbuf==NULL) {
    printf("Malloc failed to provide inbuf[%ld]\n", buf_len);
  }
  printf("inbuf has size = %ld\n", buf_len);
  unsigned char *buf_ptr = inbuf;

  printf("Generating packets for BOS\n");
  for(uint64_t i =0; i< PARAM_N; i++){ 
    packet *packet_to_process = generate_packet(params.n, params.b);
    unsigned char *packet_serialized = serialize_packet(packet_to_process);
    memcpy(buf_ptr, packet_serialized, SERIALIZED_PACKET_SIZE);
    buf_ptr+=SERIALIZED_PACKET_SIZE;
  
    free(packet_to_process);
    free(packet_serialized);  
  }  
  printf("Finished generating packets for BOS\n");

 //BucketOSort_computeParams(n, lambda, &params); 
  printf("Before BOS_Initialize call\n");

  // Time
  process_start = clock();
  BucketOSort_initialize(&params); 
  BucketOSort_sort(inbuf, inbuf, buf_len);
  process_stop = clock();
  ptime = double(process_stop-process_start)/double(CLOCKS_PER_MS);
  printf("Sort time = %f ms\n", ptime);
  free(inbuf);
  */
  return 0; 
}
