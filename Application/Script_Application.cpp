#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <openssl/err.h>
#include <math.h>
#include "Application.hpp"
#include "gcm.h"
#include <fstream>
#include <math.h>
#include <sstream>
#include <cfloat>
#include <limits>

#define NUM_ARGUMENTS_REQUIRED_1 4
#define NUM_ARGUMENTS_REQUIRED_2 7

// IV = 12 bytes, TAG = 16 bytes for AES-GCM
// So data component needs to be atleast >28 bytes large to do 0-encrypted integrity check
// We use 32 as it is the next fit for the block sizes supported by our library.
#define DATA_ENCRYPT_MIN_SIZE 40

// Global parameters that are to be supplied by the script that runs this application:
uint8_t MODE;
size_t N;
// BLOCK_SIZE > 8 (bytes), and should be a multiple of 16 bytes (For OSWAP to work) 
size_t BLOCK_SIZE;
size_t REPEAT;
size_t f_OPT=0, d_OPT=0, B_OPT=0;
bool ENCRYPTION=true;

uint64_t NUM_ZERO_BYTES;
unsigned char *zeroes;

//NOTE: Take this off if we dont care about depth
int MIN_D=1;

unsigned int countSetBits(int n) { 
  unsigned int count = 0; 
  while (n) { 
      n &= (n - 1); 
      count++; 
  } 
  return count; 
} 

void printBufferKeys(unsigned char * buf, size_t N, size_t block_size) {
  unsigned char *buf_ptr = buf;
  for(size_t i=0; i<N; i++) {
    printf("%ld, ", *((int64_t*) buf_ptr));
    buf_ptr+=block_size; 
  }
  printf("\n");
}

void printSelected(bool* selected, size_t N) {
  for(size_t i=0; i<N; i++) {
    printf("%d, ", selected[i]);
  }
  printf("\n");
}

//生成随机布尔数组
void randomizeSelectedList(bool *selected, size_t N, int randfd) {
  size_t bits_per_sizet = 8*sizeof(size_t);
  size_t bitnum = bits_per_sizet;
  size_t randword = 0;
  for(size_t i=0; i<N; i++){
    if (bitnum == bits_per_sizet) {
      read(randfd, &randword, sizeof(randword));
    }
    --bitnum;
    selected[i]=((randword >> bitnum)&1);
    if (bitnum == 0) bitnum = bits_per_sizet;
  }
}

double calculatePAPCost(size_t N, size_t B, int f, int d, size_t block_size) {
  double time;
  //Calculate time
  // time = m_1 * B + m_2 * block_size + m_3 * B * block_size + c
  // From bos_opt_mlr.py:
  // c =  0.01881 , m_1 = 9.6848e-03, m_2 = 2.6889e-06, m_3 = 7.4370e-05
  double m_1 = 9.6848e-3;
  double m_2 = 2.6889e-6;
  double m_3 = 7.4370e-5;
  double c = 0.01881;
  
  // Additional PAPs from flushBuffer execution
  // (2 * d * N for the 2 N packets passed through d sized network
  //  d * f^(d-1) * B * f from flush buffers -> d * f^d * B
  size_t N_PAP= 2*d*N + d*pow(f,d)*B;

  double PAP_cost = c + m_1 * B + m_2 * block_size + m_3 * B * block_size;
  time = N_PAP * PAP_cost;
  // Time in us (micro seconds) is returned by the formula 
  // Converting time to ms
  time/= 1000;
  return time;
}

double calculateTCCost(size_t Z_prime, int f, int d, size_t block_size) {
  double time;
  //Calculate time
  // Z'' = Z_prime log Z_prime
  size_t ZlgZ = Z_prime * ceil(log2(Z_prime));
  // time = m_1 * Z'' + m_2 * Z'' * block_size + c
  // m_1 = 6.1337-03
  // m_2 = 4.5100e-05
  // c = -2.9782
  //double m_1 = 6.1337e-3;
  //double m_2 = 4.5100e-5;
  //double c = -2.9782;


  //Removing the negative coeffient on block_size
  double m_1 =  5.9079e-3;
  double m_2 =  5.12789e-5;
  double c = -82.8964;

  size_t N_TC = pow(f,d);
  double TC_cost = (m_1 * ZlgZ) + (m_2 * ZlgZ * block_size) + c;
  //printf("For TC_cost prediction N_TC = %ld, Z_prime=%ld, each_TC_cost = %f\n", N_TC, Z_prime, TC_cost/1000);

  time = N_TC * TC_cost;
  //printf("N_TC = %ld, TC_cost = %f, total_TC_time=%f\n", N_TC, TC_cost, time/1000);
  // Time in us (micro seconds) is returned by the formula 
  // Converting time to ms
  time/=1000;
  return time;
}


double calculateRSCost(size_t N, int f, int d, size_t block_size) {
  double time;
  //Calculate time
  // N_RS = reals_per_bucket = N/b = ceil(N/(f^d))
  size_t N_buckets = pow(f,d);
  size_t N_RS = ceil((double)N / (double)N_buckets);
  
  // NlgsqN = N.log^2(N)
  size_t NlgsqN = N_RS * pow(ceil(log2(N_RS)),2);
  // time = m_1 * NlgsqN + m_2 * block_size + m_3 * NlgsqN * block_size + c
  // Values based on RS problem sizes upto 100000:
  // m_1 =  1.6558e-3
  // m_2 = -4.0110e+1
  // m_3 =  5.0619e-5
  // c =    6469.8852
  // To make it tighter, extracted values on RS problem sizes upto 50000, which gives a better fit, and is a tight upper bound to the real problem sizes we deal with
  // m_1 =  3.9360e-3
  // m_2 = -4.8680e+0
  // m_3 =  3.6198e-5
  // c =    681.317


  //Removing the negative coeffient on block_size
  // m_1 = 4.1267e-3
  // m_2 = 3.5645e-5
  // c = -997.6142

  // Reduce to <10000 sized problems on predictor
  double m_1 = 5.6801e-3;
  double m_2 = 1.9959e-5;
  double c = 145.5115;

  double RS_cost = (m_1 * NlgsqN) + (m_2 * NlgsqN * block_size) + c;
  time = N_buckets * RS_cost;
  //printf("N_RS = %ld, T1 = %f, T2 = %f, c = %f\n", N_RS, (m_1 * NlgsqN), (m_2 * NlgsqN * block_size), c);
  // Time in us (micro seconds) is returned by the formula 
  // Converting time to ms
  //printf("N_buckets = %ld, NlsqsN = %ld, RS_cost = %f, total_RS_time=%f\n", N_buckets, NlgsqN, RS_cost, time/1000);
  time/=1000;
  return time;
}
size_t calculateZ_prime(size_t N, int f, int d, int B) {
  size_t denom = pow(f, d);
  size_t Z_prime = ceil(double(2*N)/double(denom))+d*B;
  //printf("Z_prime = %ld\n", Z_prime);
  return Z_prime;
}

double calculateMean(double *input, size_t N) {
double total=0;
for(size_t i=0; i<N; i++) {
  total+=input[i];   
  } 
  return(total/N);
}

double calculateStdDev(double *input, size_t N) {
  double mean = calculateMean(input, N);
  double mean_diff = 0;
  for(size_t i=0; i<N; i++) {
    mean_diff+=(pow(input[i]-mean,2)); 
  }
  double variance = (mean_diff/N);
  return (sqrt(variance));
}

double calculateStdDevGivenMean(double *input, size_t N, double mean) {
  double mean_diff = 0;
  for(size_t i=0; i<N; i++) {
    mean_diff+=(pow(input[i]-mean, 2)); 
  }
  double variance = (mean_diff/N);
  return (sqrt(variance));
}

int testCompactionCorrectness(unsigned char *buffer, size_t N, size_t block_size, 
      unsigned char *og_buf, bool *selected){
  bool cc = true;

  unsigned char *buf_ptr = buffer;
  unsigned char *og_ptr = og_buf;
  bool *sel_ptr = selected;
  size_t cnt = 0;

  while(1){
    while(*(sel_ptr)!=1) {
      sel_ptr+=1;
      og_ptr+=block_size;
      cnt+=1;
      if(cnt==N)
        return cc;
    }

    size_t buf_key, og_key;
    buf_key = *(int64_t*) buf_ptr;
    og_key = *(int64_t*) og_ptr;
    //printf("Comparing OG: %d and buf %d (cnt = %d)\n", og_key, buf_key, cnt);
    if(buf_key != og_key) {
      //printf("Failed at index = %d\n", cnt);
      return false;
    }
    buf_ptr += block_size; 
    og_ptr += block_size;
    cnt+=1;
    sel_ptr+=1;
    if(cnt == N)
      return cc;
  }

  //Check that in buffer, for every 1 in index i of selected, og_buf[i] = buf[j] 
  return cc;
}

int testSortCorrectness(unsigned char *buffer, size_t N, size_t block_size){
  unsigned char *bfr_ptr = buffer;
  uint64_t last_key=0;
  uint64_t current_key=0;
  bool correctness = true;
  #ifdef SHOW_INPUT_KEYS
    printf("\nIn testSort: Keys after sort:\n");
  #endif
  for(size_t i=0; i<N; i++) {
    memcpy(&current_key, bfr_ptr, sizeof(uint64_t));
    #ifdef SHOW_INPUT_KEYS
      printf("%ld, ",current_key);
    #endif
    if(last_key > current_key){
      correctness = false;
      //break;
    }
    last_key = current_key;
    bfr_ptr+=block_size;
  }
  #ifdef SHOW_INPUT_KEYS
    printf("\n");
  #endif
  if(correctness == false){
    printf("Sort Correctness Failed\n");
    return -1;
    #ifdef SHOW_INPUT_KEYS
      printf("\n\n");
    #endif
  }
  #ifdef SHOW_INPUT_KEYS
    printf("\n\n");
  #endif
  
  return 0;
}

int testShuffleCorrectness(unsigned char *buffer, size_t N, size_t block_size, int *key_counts){
  unsigned char *bfr_ptr = buffer;
  int *returned_key_counts = new int[N]{};
  uint64_t current_key;   
  
  #ifdef SHOW_INPUT_KEYS
    printf("In testShuffle: Keys after shuffle:\n");
  #endif
  for(size_t i=0; i<N; i++) {
    memcpy(&current_key, bfr_ptr, sizeof(uint64_t));
    #ifdef SHOW_INPUT_KEYS
      printf("%ld, ",current_key);
    #endif
    if(current_key>=N){
      printf("Got a >= N key!!!!!!!!\n");
    } else {
      returned_key_counts[current_key]+=1;
    }
    bfr_ptr+=block_size;
  }
  #ifdef SHOW_INPUT_KEYS
    printf("\n\n\n");
  #endif

  int memcmp_ret = memcmp(key_counts, returned_key_counts, N);
  //delete [] returned_key_counts;

  if(memcmp_ret){
    printf("Shuffle Correctness Failed\n");
    return 1;
  }
 
  return 0;  
}

void printBOSParams(bos_params *params) {
  printf("BOS Params: \n");
  printf("N = %ld\n", params->n);
  printf("f = %d\n", params->f);
  printf("b = %ld\n", params->b);
  printf("d = %d\n", params->d);
  printf("Z = %ld\n", params->Z);
  printf("B = %ld\n", params->B);
  printf("block_size = %ld\n\n", params->block_size);
}

void parseCommandLineArguments(int argc, char *argv[]){
  if(!(argc==(NUM_ARGUMENTS_REQUIRED_1 + 1) or argc==(NUM_ARGUMENTS_REQUIRED_2 + 1)))  {
    printf( "Did NOT receive the right number of command line arguments.\n"
      "Style (1): ./script_application <MODE> <N> <BLOCK_SIZE> <REPEAT>\n"
      "Style (2): ./script_application <3> <N> <BLOCK_SIZE> <REPEAT> <f_opt> <d_opt> <B_opt>\n\n"
      "Style 2 is for BORP Modes (3 and 13), where <f_opt>,<d_opt>,<B_opt>\n"
      "are the BRN parameters computed by bos_optimizer.py for a given problem\n"
      "(N, BLOCK_SIZE) and the optimization goal.\n\n"

      "NOTE:\n   i) BLOCK_SIZE should be 8, or a multiple of 16! \n"
      "  ii) For Sorts, the first 8 bytes of a block are used as its label, \n"
      "      except for MODE 17 which uses just the first 4 bytes as the label.\n\n");

    printf("MODE LIST:\n"
            "  Oblivious Shuffling Algorithms (1 to 9)\n"
            "    (1) RecursiveShuffle (with Goodrich Compaction)\n"
            "    (2) ORShuffle\n"
            "    (3) BORPStream (Bucket Oblivious Random Permutation Stream) (lambda = -80)\n"
            "        (NOTE: BORP Modes (3 and 13) should be run through run_experiments.py\n"
            "        or/i.e. by running bos_optimizer.py first and then using the script in\n"
            "        the second use style.)\n"
            "    (4) OddEvenMergeSort Shuffle\n"
            "    (6) ButterflyNetwork Compaction (新增)\n"  
            "    (5) BitonicShuffle\n"
            "    (7) WaksShuffle\n"
            "    (8) Bucket Oblivious Random Permutation (with TC)\n"
            "    (9) Nassimi-Sahni Waksman Shuffle\n\n"

            "  Oblivious Sorting Algorithms (10 to 17)\n"
            "    (10) Nassimi-Sahni Waksman Sort\n"
            "    (11) Oblivious Sort (Sorting Network - OddEvenMergeSort)\n"
            "    (12) Oblivious Sort (Sorting Network - BitonicSort)\n"
            "    (13) Oblivious Sort (BORPStream + Quicksort, lambda = -80)\n"
            "        (NOTE: BORP Modes (3 and 13) should be run through run_experiments.py\n"
            "        or/i.e. by running bos_optimizer.py first and then using the script in\n"
            "        the second use style.)\n"
            "    (14) Oblivious Sort (ORShuffle + Quicksort)\n"
            "    (15) WaksShuffle + Quicksort\n"
            "    (17) WaksSort\n\n"

            "  Oblivious Compaction Algorithms (21 & 22)\n"
            "    (21) ORCompact\n"
            "    (22) Goodrich Compaction\n");
            //"61 = RecursiveShuffle_M2\n"
            //"91 = MeasureOSWAPBuffer\n");
   exit(0);
  }

  MODE = atoi(argv[1]);
  N = atoi(argv[2]);
  BLOCK_SIZE = atoi(argv[3]);
  REPEAT = atoi(argv[4]);
  // To ignore the first iteration, we perform the experiment REPEAT + 1 times
  REPEAT = REPEAT + 1;
  if(argc==(NUM_ARGUMENTS_REQUIRED_2+1)) {
    f_OPT = atoi(argv[5]);
    d_OPT = atoi(argv[6]);
    B_OPT = atoi(argv[7]);
  }
}

// Return real times in microseconds
uint64_t rtclock() {
    static time_t secstart = 0;
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    if (secstart == 0) secstart = tp.tv_sec;
    return (tp.tv_sec-secstart) * 1000000 + tp.tv_nsec / 1000;
}

int main(int argc, char *argv[]){
  clock_t process_start, process_stop, outer_encryption_start, outer_encryption_stop, 
    outer_decryption_start, outer_decryption_stop;
  clock_t process_post_init;
  double itime, ptime;
  double pure_ptime;  
  size_t nthreads = 1;
  char *threadenv = getenv("ENCLAVE_THREADS");
  if (threadenv) {
    int n = atoi(threadenv);
    if (n >= 1) {
        nthreads = n;
    }
  }

  bool verbose_phases = !!getenv("VERBOSE_PHASES");
  uint64_t phase_start = rtclock(), phase_end;
  double phase_time;

  // Load experiment parameters
  parseCommandLineArguments(argc, argv);

  phase_end = rtclock();
  phase_time = (double)(phase_end - phase_start)/1000.0;
  if (verbose_phases) {
      printf("parseCommandLineArguments phase: %f\n", phase_time);
  }
  phase_start = phase_end;

  double ptime_array[REPEAT] = {};
  double pure_ptime_array[REPEAT] = {};
  // For Waksman, num_oswaps returns the OSWAPS performed while applying permutation
  // and num_owaps_cb (below) returns the number of OSWAPs performed during the control
  // bit setting algorithm
  size_t num_oswaps[REPEAT] = {};
  double lltime[REPEAT]= {};
  size_t num_oswaps_gp[REPEAT] = {};
  size_t num_oswaps_cb[REPEAT] = {};
  size_t num_oswaps_ap[REPEAT] = {};
  double waksman_rp_time[REPEAT]={};
  double waksman_cb_time[REPEAT]={};
  double waksman_ap_time[REPEAT]={};
  double qsort_time[REPEAT]={};

  double butterfly_cb_time[REPEAT]={};
  double butterfly_ap_time[REPEAT]={};


  // Initialize libcrypto
  OpenSSL_add_all_algorithms();
  ERR_load_crypto_strings();

  phase_end = rtclock();
  phase_time = (double)(phase_end - phase_start)/1000.0;
  if (verbose_phases) {
      printf("initalize libcrypto phase: %f\n", phase_time);
  }
  phase_start = phase_end;

  int randfd = open("/dev/urandom", O_RDONLY);
  if (randfd < 0) {
    throw std::runtime_error("Cannot open /dev/urandom");
  }

  // Initialize enclave
  OLib_initialize();

  phase_end = rtclock();
  phase_time = (double)(phase_end - phase_start)/1000.0;
  if (verbose_phases) {
      printf("OLib_initialize phase: %f\n", phase_time);
  }
  phase_start = phase_end;

  // Load test AES keys into the enclave
  unsigned char inkey[16];
  unsigned char outkey[16];
  unsigned char datakey[16];
  read(randfd, inkey, 16);
  read(randfd, outkey, 16);
  read(randfd, datakey, 16);
  Enclave_loadTestKeys(inkey, outkey);

  const size_t ENC_BLOCK_SIZE = 12 + BLOCK_SIZE + 16; 

  phase_end = rtclock();
  phase_time = (double)(phase_end - phase_start)/1000.0;
  if (verbose_phases) {
      printf("loadtestkeys phase: %f\n", phase_time);
  }
  phase_start = phase_end;

  // Create buffer of items to shuffle
  
  if(MODE==21 or MODE==22) {
    ENCRYPTION=false; 
  }

  unsigned char *buf;
  // For testing compact algs we need the original buffer and selected list at the end
  unsigned char *og_buf;
  size_t buflen, decbuflen;

  if(ENCRYPTION){
    buflen = N * ENC_BLOCK_SIZE;
  } else {
    buflen = N * BLOCK_SIZE;
  }
  buf = new unsigned char[buflen];

  if (buf == NULL) {
    printf("Allocating buffer memories in script Application failed!\n");
  }

  unsigned char *bufend = buf + buflen;
  ptime = 0;

  phase_end = rtclock();
  phase_time = (double)(phase_end - phase_start)/1000.0;
  if (verbose_phases) {
      printf("selected_list phase: %f\n", phase_time);
  }
  phase_start = phase_end;


  bool *selected_list = NULL;
  if(MODE==21||MODE==22)
    selected_list = new bool[N]{};
  int *key_counts;

  for (size_t r=0; r<REPEAT; r++) {

    if(MODE>=1 && MODE <=9) {
      key_counts = new int[N]{};
    }
    
    if(MODE==21 || MODE ==22) {
      randomizeSelectedList(selected_list, N, randfd);
    }
    
    size_t inc_ctr=0;
    if(ENCRYPTION) {
      unsigned char iv[12];
      read(randfd, iv, 12);
      for (unsigned char *enc_block_ptr = buf; enc_block_ptr < bufend; enc_block_ptr += ENC_BLOCK_SIZE) {
        unsigned char block[BLOCK_SIZE] = {}; // Initializes to zero


        uint64_t rnd;
        #ifdef RANDOMIZE_INPUTS
          read(randfd, (unsigned char*) &rnd, sizeof(rnd));
          //For easier visual debugging:
          rnd = rnd % N;
        #else
          rnd = inc_ctr++;
        #endif

        if(MODE>=1 && MODE <=9) {
          key_counts[rnd]+=1;
        }
        
        #ifdef SHOW_INPUT_KEYS
          printf("%ld, ",rnd);
        #endif

        if(MODE==7 || MODE==17 || MODE ==15) {
          uint32_t rnd32 = rnd;
          memcpy(block, (unsigned char*) &rnd32, sizeof(rnd32));
        }
        else{
          memcpy(block, (unsigned char*) &rnd, sizeof(rnd));
        }

   
        if(BLOCK_SIZE >= DATA_ENCRYPT_MIN_SIZE) {
          // NUM_ZERO_BYTES = the number of zero bytes that get encrypted
          // i.e. BLOCK_SIZE - 8 (key bytes) - 12 (IV bytes) - 16 (TAG bytes)
          NUM_ZERO_BYTES = BLOCK_SIZE-8-12-16;
          zeroes = new unsigned char[NUM_ZERO_BYTES]();
        
          /*
          printf("zeroes = ");
          for(int p=0; p<NUM_ZERO_BYTES; p++) {
            printf("%x,", zeroes[p]);
          }
          printf("\n");
          */

          (*((uint64_t*)iv))++;
          memmove(block+8, iv, 12);

          if ((NUM_ZERO_BYTES) != gcm_encrypt(zeroes, NUM_ZERO_BYTES, NULL, 0, datakey,
                  block+8, 12, block+8+12, block+8+12+NUM_ZERO_BYTES)) {
              printf("Encryption failed\n");
              break;
          }/* else{
            printf("Encrypted_data_block = ");
            for(int p=0; p<BLOCK_SIZE-8; p++) {
              printf("%x", block[8+p]);
            }
            printf("\n");
          }*/  
        }

        // Encrypt the chunk to the enclave      
        outer_encryption_start = clock();
        (*((uint64_t*)iv))++;
        memmove(enc_block_ptr, iv, 12);
        if (BLOCK_SIZE != gcm_encrypt(block, BLOCK_SIZE, NULL, 0, inkey,
                enc_block_ptr, 12, enc_block_ptr+12, enc_block_ptr + 12 + BLOCK_SIZE)) {
            printf("Encryption failed\n");
            break;
        }
        outer_encryption_stop = clock();
        itime = double(outer_encryption_stop-outer_encryption_start)/double(CLOCKS_PER_MS);
        ptime+=itime;
      }
      #ifdef SHOW_INPUT_KEYS
        printf("\n");
      #endif
      // Print outer_encryption_time
      //printf("%f\n", ptime);
      // NOTE: Mute the below line. (Meant for visual key debugging purposes)
      //printf("\n\n");

      phase_end = rtclock();
      phase_time = (double)(phase_end - phase_start)/1000.0;
      if (verbose_phases) {
          printf("preparation phase: %f\n", phase_time);
      }

    } else {
      uint64_t rnd;
      for (unsigned char *block_ptr = buf; block_ptr < bufend; block_ptr += BLOCK_SIZE) {
        unsigned char block[BLOCK_SIZE] = {}; // Initializes to zero

        //read(randfd, (unsigned char*) &rnd, sizeof(rnd));
        //For easier visual debugging:
        //rnd = rnd % N;
        rnd = inc_ctr++;

        if(MODE>=1 && MODE <=9) {
          key_counts[rnd]+=1;
        }
        
        //printf("%ld, ",rnd);
        memcpy(block, (unsigned char*) &rnd, sizeof(rnd));
        memcpy(block_ptr, block, BLOCK_SIZE);
      }
      phase_end = rtclock();
      phase_time = (double)(phase_end - phase_start)/1000.0;
      if (verbose_phases) {
          printf("preparation phase: %f\n", phase_time);
      }
      size_t buflen = N * BLOCK_SIZE;
      og_buf = new unsigned char[buflen];
      memcpy(og_buf, buf, buflen); 
    } 
    phase_start = phase_end;
   
  
    process_start = clock();
    // Select operation based on MODE:
  
    // For the BOS modes we need a BRN_params structure that gets populated by BOS_optimize
    BRN_params params;
    bos_params bosparams;
    enc_ret ret;

    switch(MODE){
      case 1:
        pure_ptime = DecryptAndShuffleM1(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;  

      case 2:
        pure_ptime = DecryptAndShuffleM2(buf, N, ENC_BLOCK_SIZE, nthreads, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;

      case 3: 
        //NOTE: BORP Streaming expects the f_OPT, d_OPT, and B_opt params to be supplied from the optimizer (bos_optimizer.py)
        bosparams.n = N;
        bosparams.n_1 = N;
        bosparams.block_size = BLOCK_SIZE;
        bosparams.evict_start = 0;
        bosparams.f = f_OPT;
        bosparams.d = d_OPT;
        bosparams.B = B_OPT;

        bosparams.b = pow(f_OPT, d_OPT);
        bosparams.Z = ceil(2 * (double)N/(double(bosparams.b)));

        //printBOSParams(&bosparams);
        pure_ptime = DecryptAndBORP(buf, N, ENC_BLOCK_SIZE, &bosparams, nthreads, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        lltime[r]=ret.last_layer_time;
        break;

      case 4:
        pure_ptime = DecryptAndOddEvenMergeSortShuffle(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;
      
      case 5:
        pure_ptime = DecryptAndBitonicSortShuffle(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;

      case 6:
        DecryptAndOblivButterflyCompact(buf, (uint32_t) N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r] = ret.ptime;
        num_oswaps[r] = ret.OSWAP_count;
        num_oswaps_cb[r] = ret.OSWAP_cb;    // 控制位计算阶段的交换次数
        num_oswaps_ap[r] = ret.OSWAP_ap;    // 应用阶段的交换次数
        butterfly_cb_time[r] = ret.control_bits_time;  // 控制位计算时间
        butterfly_ap_time[r] = ret.apply_perm_time;    // 应用阶段时间
        break;

      case 7:
        DecryptAndOblivWaksmanShuffle(buf, (uint32_t) N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=ret.ptime;
        num_oswaps[r]=ret.OSWAP_ap;
        num_oswaps_gp[r]= ret.OSWAP_gp;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        waksman_rp_time[r] = ret.gen_perm_time;
        waksman_cb_time[r] = ret.control_bits_time;
        waksman_ap_time[r] = ret.apply_perm_time;
        break;

      case 8:
        pure_ptime = DecryptAndBORP_TC(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;

      case 9:
        DecryptAndDJBWaksmanShuffle(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=ret.ptime;
        num_oswaps[r]=ret.OSWAP_ap;
        num_oswaps_gp[r]= ret.OSWAP_gp;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        waksman_rp_time[r] = ret.gen_perm_time;
        waksman_cb_time[r] = ret.control_bits_time;
        waksman_ap_time[r] = ret.apply_perm_time;
        break;

      case 10:
        DecryptAndDJBWaksmanSort(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=ret.ptime;
        num_oswaps[r]=ret.OSWAP_ap;
        num_oswaps_gp[r]= ret.OSWAP_gp;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        waksman_rp_time[r] = ret.gen_perm_time;
        waksman_cb_time[r] = ret.control_bits_time;
        waksman_ap_time[r] = ret.apply_perm_time;
        break;

      case 11:
        pure_ptime = DecryptAndOddEvenMergeSort(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;
      
      case 12:
        pure_ptime = DecryptAndBitonicSort(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;
      
      case 13: 
        //NOTE: BORP Streaming expects the f_OPT, d_OPT, and B_opt params to be supplied from the optimizer (bos_optimizer.py)
        bosparams.n = N;
        bosparams.n_1 = N;
        bosparams.block_size = BLOCK_SIZE;
        bosparams.evict_start = 0;
        bosparams.f = f_OPT;
        bosparams.d = d_OPT;
        bosparams.B = B_OPT;

        bosparams.b = pow(f_OPT, d_OPT);
        bosparams.Z = ceil(2 * (double)N/(double(bosparams.b)));

        //printBOSParams(&bosparams);
        pure_ptime = DecryptAndBOS(buf, N, ENC_BLOCK_SIZE, &bosparams, nthreads, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        lltime[r]=ret.last_layer_time;
        qsort_time[r]=ret.qsort_time;
        break;
      
      case 14:
        pure_ptime = DecryptAndOSortwithRSandNS(buf, N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        qsort_time[r]=ret.qsort_time;
        break;

      case 15:
        //DecryptAndObliviousWakshmanShuffle + qSort
        DecryptAndOWSS(buf, (uint32_t) N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=ret.ptime;
        num_oswaps[r]=ret.OSWAP_ap;
        num_oswaps_gp[r]= ret.OSWAP_gp;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        waksman_rp_time[r] = ret.gen_perm_time;
        waksman_cb_time[r] = ret.control_bits_time;
        waksman_ap_time[r] = ret.apply_perm_time;
        qsort_time[r] = ret.qsort_time;
        break;

      case 17:
        DecryptAndOblivWaksmanSort(buf, (uint32_t) N, ENC_BLOCK_SIZE, buf, &ret);
        pure_ptime_array[r]=ret.ptime;
        num_oswaps[r]=ret.OSWAP_ap;
        num_oswaps_gp[r]= ret.OSWAP_gp;
        num_oswaps_cb[r]=ret.OSWAP_cb;
        num_oswaps_ap[r]=ret.OSWAP_ap;
        waksman_rp_time[r] = ret.gen_perm_time;
        waksman_cb_time[r] = ret.control_bits_time;
        waksman_ap_time[r] = ret.apply_perm_time;
        qsort_time[r] = ret.qsort_time;
        break;

      case 21:
        //printBufferKeys(buf, N, BLOCK_SIZE);
        //printSelected(selected_list, N);
        printf("Before TTC\n");
        pure_ptime = testTightCompaction(buf, N, BLOCK_SIZE, nthreads, selected_list, &ret);
        printf("After TTC\n");
        //printBufferKeys(buf, N, BLOCK_SIZE);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break;


      case 22:
        pure_ptime = testOPTightCompaction(buf, N, BLOCK_SIZE, selected_list, &ret);
        pure_ptime_array[r]=pure_ptime;
        num_oswaps[r]=ret.OSWAP_count;
        break; 

      case 61:
        // NOTE we don't care about correctness here, we are sending oversized encrypted packets
        // and treating it as plaintext packets, just to get the raw timing of RS without
        // encryption and decryption overheads.
        pure_ptime = RecursiveShuffle_M2(buf, N, BLOCK_SIZE);
        pure_ptime_array[r] = pure_ptime;
        break;

      case 91:
        // It performs OSWAPs on the N blocks of buf, with a local "result" buf,
        // the swap flag is set by the first bit of each buf
        // NOTE: We don't care that the buffer is encrypted, and has N blocks of ENC_BLOCK_SIZE
        // This function is only designed to test timings of OSWAP_BUFFER calls (like cmp vs test)
        pure_ptime = MeasureOSWAPBuffer(buf, N, BLOCK_SIZE);
        pure_ptime_array[r]=pure_ptime;
        break;

      case 92:
        // It performs OSWAPs on the N blocks of buf in a pairwise fashion
        //
        break;
    }
    process_stop = clock();

    ptime = double(process_stop-process_start)/double(CLOCKS_PER_MS);
    ptime_array[r]=ptime; 

    phase_end = rtclock();
    phase_time = (double)(phase_end - phase_start)/1000.0;
    if (verbose_phases) {
        printf("processing phase: %f\n", phase_time);
    }
    phase_start = phase_end;

    //TODO: Remove 13 from exception list!
    bool dec_fail_flag = false;
    if(ENCRYPTION) {
      outer_decryption_start = clock();
      int numfailed = 0;
      int cnum = 0;
       
      
      unsigned char *decrypted_result_buf_ptr = buf; 
      for (unsigned char *enc_block_ptr = buf; enc_block_ptr < bufend; enc_block_ptr += ENC_BLOCK_SIZE) {
        unsigned char block[BLOCK_SIZE];
        ++cnum;
        if (BLOCK_SIZE != gcm_decrypt(enc_block_ptr+12, BLOCK_SIZE, NULL, 0, enc_block_ptr+12+BLOCK_SIZE, outkey, enc_block_ptr, 12, block)) {
            printf("Outer Decryption failed %d/%d\n", ++numfailed, cnum);
            dec_fail_flag = true;
            break;
        }


        if(BLOCK_SIZE >= DATA_ENCRYPT_MIN_SIZE) {
          //Check correctness of data_payload
          unsigned char should_be_zeroes[NUM_ZERO_BYTES]={};
          /*
          printf("Output Encrypted_data_block = ");
          for(int p=0; p<BLOCK_SIZE-8; p++) {
            printf("%x", block[8+p]);
          }
          printf("\n");
          */
          int returned_pt_bytes = gcm_decrypt(block+8+12, NUM_ZERO_BYTES, NULL, 0,
                  block+8+12+NUM_ZERO_BYTES, datakey, block+8, 12, should_be_zeroes); 
          if (NUM_ZERO_BYTES != returned_pt_bytes) {
            printf("Data block decryption failed, returned_pt_bytes = %d\n", returned_pt_bytes);  
            dec_fail_flag = true;
            break;
          }/*
          else{
            if(memcmp(should_be_zeroes, zeroes, NUM_ZERO_BYTES)) {
              
              printf("should_be_zeroes = ");
              for(int p=0; p<NUM_ZERO_BYTES; p++) {
                printf("%x,", should_be_zeroes[p]);
              }
              printf("\n");
            }
          }         */
        }

        memcpy(decrypted_result_buf_ptr, block, BLOCK_SIZE);
        decrypted_result_buf_ptr+=BLOCK_SIZE;
      }
      outer_decryption_stop = clock();
      ptime = double(outer_decryption_stop - outer_decryption_start)/double(CLOCKS_PER_MS);
      if(dec_fail_flag){
        exit(0);
      }
    }      
       
    switch(MODE){
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 7:
      case 8: 
      case 9: {
                testShuffleCorrectness(buf, N, BLOCK_SIZE, key_counts);
                break;
              }
      case 10: 
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:{
                testSortCorrectness(buf, N, BLOCK_SIZE);
                break;
              }
      case 21:
      case 22:{

                printf("Before TC Correctness check\n");
                /*
                bool tcc_flag = testCompactionCorrectness(buf, N, BLOCK_SIZE, og_buf, selected_list);
                if(!tcc_flag){
                  printf("Order-preserving compaction INCORRECT\n");
                }
                */
                printf("After TC Correctness check\n");
                break;
              }
    } 
    phase_end = rtclock();
    phase_time = (double)(phase_end - phase_start)/1000.0;
    if (verbose_phases) {
        printf("check phase: %f\n", phase_time);
    }
    phase_start = phase_end;
  }

  //NOTE: The +1 and -1 are to remove the additional timing variance that stems from the 
  //      first run of execution always taking longer to execute due to instruction page
  //      cache warming.
  double mean = calculateMean(ptime_array+1, REPEAT-1);
  double stddev = calculateStdDevGivenMean(ptime_array+1, REPEAT-1, mean);
  double llmean = calculateMean(lltime+1, REPEAT-1);
  double llstddev = calculateStdDevGivenMean(lltime+1, REPEAT-1, llmean);

  double mean_pure = calculateMean(pure_ptime_array+1, REPEAT-1);
  double stddev_pure = calculateStdDevGivenMean(pure_ptime_array+1, REPEAT-1, mean_pure);
  printf("%f\n", mean);
  printf("%f\n", stddev);
  printf("%f\n", mean_pure);
  printf("%f\n", stddev_pure);
  printf("%ld\n", num_oswaps[0]);

  if(MODE==7 || MODE==9 || MODE==10 || MODE==15 || MODE ==17){
    double rp_mean = calculateMean(waksman_rp_time+1, REPEAT-1);
    double cb_mean = calculateMean(waksman_cb_time+1, REPEAT-1);
    double ap_mean = calculateMean(waksman_ap_time+1, REPEAT-1);
    double rp_stddev = calculateStdDevGivenMean(waksman_rp_time + 1, REPEAT-1, rp_mean);
    double cb_stddev = calculateStdDevGivenMean(waksman_cb_time + 1, REPEAT-1, cb_mean);
    double ap_stddev = calculateStdDevGivenMean(waksman_ap_time + 1, REPEAT-1, ap_mean);
    printf("%f\n", rp_mean);
    printf("%f\n", rp_stddev);
    printf("%f\n", cb_mean);
    printf("%f\n", cb_stddev);
    printf("%f\n", ap_mean);
    printf("%f\n", ap_stddev);
    printf("%ld\n", num_oswaps_gp[0]);
    printf("%ld\n", num_oswaps_cb[0]);
    printf("%ld\n", num_oswaps_ap[0]);
  } else{
    printf("%f\n", llmean);
    printf("%f\n", llstddev);
  }

  if(MODE==3 || MODE ==13) {
    printf("%ld\n", num_oswaps_cb[0]);
    printf("%ld\n", num_oswaps_ap[0]);
  }

  if( MODE==13 || MODE==14 || MODE ==15){
    double qsort_mean = calculateMean(qsort_time+1, REPEAT-1);
    double qsort_stddev = calculateStdDevGivenMean(qsort_time+1, REPEAT-1, qsort_mean);
    printf("%f\n", qsort_mean);
    printf("%f\n", qsort_stddev);
  }

  if (MODE == 6) {
    double cb_mean = calculateMean(butterfly_cb_time+1, REPEAT-1);
    double ap_mean = calculateMean(butterfly_ap_time+1, REPEAT-1);
    double cb_std  = calculateStdDevGivenMean(butterfly_cb_time+1, REPEAT-1, cb_mean);
    double ap_std  = calculateStdDevGivenMean(butterfly_ap_time+1, REPEAT-1, ap_mean);

    printf("%f\n", cb_mean);
    printf("%f\n", cb_std);
    printf("%f\n", ap_mean);
    printf("%f\n", ap_std);
    printf("%ld\n", num_oswaps_cb[0]);
    printf("%ld\n", num_oswaps_ap[0]);
  }

  
  close(randfd);

  if(MODE>=1 && MODE<=9){
    //delete []key_counts;
  }

  //delete[] buf;
  if(BLOCK_SIZE >= DATA_ENCRYPT_MIN_SIZE)
    delete []zeroes;

  if(MODE==21||MODE==22){
    delete []selected_list;
    delete []og_buf;
  }
  return 0;
}
