#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <openssl/err.h>
#include <stdint.h>
#include "Application.hpp"
#include "gcm.h"
#include <fstream>
#include <math.h>
#include <sstream>
#include <cfloat>
#include <limits>

#define NUM_ARGUMENTS_REQUIRED 3

size_t N;
// BLOCK_SIZE > 8 (bytes), and should be a multiple of 16 bytes (For OSWAP to work) 
size_t BLOCK_SIZE;
size_t REPEAT;

int MIN_D = 1;

unsigned int countSetBits(int n) { 
  unsigned int count = 0; 
  while (n) { 
      n &= (n - 1); 
      count++; 
  } 
  return count; 
} 

void parseCommandLineArguments(int argc, char *argv[]){
  if(argc!=(NUM_ARGUMENTS_REQUIRED + 1))  {
    printf("Did not receive the right number of command line arguments.\n");
    printf("Expected Use: ./bos_application <N> <BLOCK_SIZE> <REPEAT>\n");
    printf("NOTE: BLOCK_SIZE should be 8, or a multiple of 16! \n");
   exit(0);
  }

  N = atoi(argv[1]);
  BLOCK_SIZE = atoi(argv[2]);
  REPEAT = atoi(argv[3]);
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
  size_t N_PAP=2*d*N + (d * (d + 1))*(pow(f,d)*B)/2;

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


// Assumes OLib_initialize call is already done in the main()
double executeBOS(BRN_params *params) {
  clock_t process_start, process_stop;
  double total_time=0, ptime;

  int randfd = open("/dev/urandom", O_RDONLY);
  if (randfd < 0) {
    throw std::runtime_error("Cannot open /dev/urandom");
  }

  // Load test AES keys into the enclave
  unsigned char inkey[16];
  unsigned char outkey[16];
  read(randfd, inkey, 16);
  read(randfd, outkey, 16);

  // An AES key for constructing and verifying the _plaintext_ of the blocks
  // so that we can ensure that the blocks come out unaltered
  unsigned char datakey[16];
  read(randfd, datakey, 16);
  
  const size_t enc_block_size = 12 + BLOCK_SIZE+ 16;

  unsigned char *inoutbuf = (unsigned char*) malloc(N * enc_block_size);
  unsigned char *outbuf = (unsigned char*) malloc(N * enc_block_size);
  unsigned char *bufend = inoutbuf + N * enc_block_size;
  unsigned char *bufout_end = outbuf + N * enc_block_size;
  size_t buf_len = N * enc_block_size;

  for (unsigned char *enc_block_ptr = inoutbuf; enc_block_ptr < bufend; enc_block_ptr+=enc_block_size) {
    unsigned char block[BLOCK_SIZE] = {}; // Initializes to zero

    uint64_t rnd;
    read(randfd, (unsigned char*) &rnd, sizeof(rnd));
    //For easier visual debugging:
    rnd = rnd % N;
    memcpy(block, (unsigned char*) &rnd, sizeof(rnd));    
    //printf("%ld, ", rnd);

    // Encrypt the chunk to the enclave      
    read(randfd, enc_block_ptr, 12); // Choose a random IV
    if (BLOCK_SIZE != gcm_encrypt(block, BLOCK_SIZE, NULL, 0, inkey,
            enc_block_ptr, 12, enc_block_ptr+12, enc_block_ptr + 12 + BLOCK_SIZE)) {
        printf("Encryption failed\n");
        break;
    }
  } 
  //printf("\n");
  close(randfd); 

  Enclave_loadTestKeys(inkey, outkey);

  bos_params bosparams;
  bosparams.f = params->f;
  bosparams.n = N;
  bosparams.n_1 = N;
  bosparams.b = pow(params->f, params->d);
  bosparams.Z = ceil(2 * (double)N/(double(bosparams.b)));
  bosparams.d = params->d;
  bosparams.B = params->B;
  bosparams.block_size = BLOCK_SIZE;
  bosparams.evict_start = 0;
  enc_ret ret;

  //BOS_initialize(&bosparams);
  size_t nthreads = 1;
  double pure_ptime = DecryptAndBOS(inoutbuf, N, enc_block_size, &bosparams, nthreads, outbuf, &ret);
  printf("Total_ptime for BOS = %f\n", pure_ptime);

  free(inoutbuf); 
  free(outbuf);
  return pure_ptime;
}


// lambda = -40, -60, or -80
void optimize(size_t N, size_t block_size, int lambda, BRN_params *opt_params){

  //File name = fdBgrdi + str(lambda)
  std::string file_base="fdBgrid";
  std::string file_suffix=std::to_string(lambda);
  std::string file_name = file_base+file_suffix;
  printf("File_name = %s\n", file_name.c_str());
  std::fstream fs;
  std::string line;

  
  fs.open(file_name.c_str(), std::fstream::in);

  //Check num lines in file
  // -1 to account for the extra count of last line from the \n
  int num_opt_points=-1;
  int min_t_index=0;
  double PAP_cost, TC_cost, RS_cost;  
  double min_t=-1;

  std::string s;
  while(!fs.eof()) {
    getline(fs, s);
    num_opt_points++; 
  }

  printf("Num of opt points = %d \n", num_opt_points); 
  fs.close();
  
  // create f[], d[], B[], t[i]
  int f[num_opt_points], d[num_opt_points], B[num_opt_points], t[num_opt_points]={};

  fs.open(file_name.c_str(), std::fstream::in);
  size_t tkn_ctr=0;
  char delim = ',';
  
  while(std::getline(fs, line)) { 
    std::string token;
    
    std::stringstream ss(line);
    int index=0;    

    while(std::getline(ss, token, delim)){
      index = tkn_ctr/4;

      if(tkn_ctr%4==0){
        f[index]=std::stoi(token, nullptr);
        //printf("f=%d, ",f[index]);
      } else if(tkn_ctr%4==1){
        d[index]=std::stoi(token, nullptr);
        //printf("d=%d, ",d[index]);
      }
      else if(tkn_ctr%4==2){
        B[index]=std::stoi(token, nullptr);
        //printf("B=%d, ",B[index]);
      }
      else if(tkn_ctr%4==3){
        //The actual FP for this configuration is populated in token in scientific notation!
        //printf("\n");
      }
      tkn_ctr+=1;
    }
 
  }   
  fs.close();

  int current_d = 0; 
  double executed_min_time = DBL_MAX;
  size_t executed_min_index = -1; 
  double feasible_min_time = DBL_MAX;
  size_t feasible_min_index = 0;

  for(int i=0; i<num_opt_points; i++) {
    size_t Z_prime = calculateZ_prime(N, f[i], d[i], B[i]);
    printf("\n\nf = %d, d = %d, B = %d:\n", f[i], d[i], B[i]);
    PAP_cost = calculatePAPCost(N, B[i], f[i], d[i], block_size);
    printf("Prediction Timing: PAP_cost = %f, ", PAP_cost);
    TC_cost = calculateTCCost(Z_prime, f[i], d[i], block_size);
    printf("TC_cost = %f, ", TC_cost);
    RS_cost = calculateRSCost(N, f[i], d[i], block_size);
    printf("RS_cost = %f, ", RS_cost);
    printf("RF_cost = %f\n", TC_cost + RS_cost);


    double total_time = PAP_cost+TC_cost+RS_cost;
    printf("Total_time = %f\n", total_time);
    t[i] = total_time; 

    if(min_t <0) {
      if(d[i]>=MIN_D)
        min_t = total_time;
        min_t_index = i;
    } else if(total_time < min_t && d[i]>= MIN_D) {

      current_d = d[i];
      min_t = total_time;
      min_t_index = i;
      printf("PAP_cost = %f, TC_cost = %f, RS_cost =%f, total_time = %f\n", PAP_cost, TC_cost, RS_cost, total_time);
    }

    if(countSetBits(f[i])==1) {
      if(total_time < feasible_min_time){
        feasible_min_index = i;
        feasible_min_time = total_time;
      }
    } 
    
    //If f is a power of 2, actually execute the BOS for this configuration
    // and compare predicted value vs actual time taken.
    if(countSetBits(f[i])==1 && f[i]<16 && d[i]<=3 && B[i]<=30) {
      BRN_params params;
      params.f = f[i];
      params.d = d[i];
      params.B = B[i];

      double exec_time = executeBOS(&params);
      
      if(executed_min_time == DBL_MAX){
        executed_min_time = exec_time;
        executed_min_index = i;
      }
      else if(exec_time < executed_min_time)  {
        executed_min_time = exec_time;
        executed_min_index = i;
      }
    } 
  }

  
  //printf("min_t_index = %d\n", min_t_index);
  // Set opt_params based on min_t_index;
  printf("\n\n\nMinimum found by Optimizer, f = %d, B = %d, d = %d\n", f[min_t_index], B[min_t_index], d[min_t_index]);
  printf("Min time = %f ms\n", min_t);
  //printf("Minimum cost: PAP = %f, TC=%f, total=%f\n",);

  printf("\n\nFEASIBLE (f = power of 2) Minimum found by Optimizer, f = %d, B = %d, d = %d\n", f[feasible_min_index], B[feasible_min_index], d[feasible_min_index]);
  printf("Min time = %f ms\n", feasible_min_time);
  opt_params->f=f[feasible_min_index];
  opt_params->d=d[feasible_min_index];
  opt_params->B=B[feasible_min_index]; 

  
  printf("\n\nMinimum found by Executions, f = %d, B = %d, d = %d\n", f[executed_min_index], B[executed_min_index], d[executed_min_index]);
  printf("Min time = %f ms\n", executed_min_time); 
} 



int main(int argc, char *argv[]){

  parseCommandLineArguments(argc, argv);

  // Initialize libcrypto
  OpenSSL_add_all_algorithms();
  ERR_load_crypto_strings();

  clock_t process_start, process_stop;
  double total_time=0, ptime;

  int randfd = open("/dev/urandom", O_RDONLY);
  if (randfd < 0) {
    throw std::runtime_error("Cannot open /dev/urandom");
  }

  printf("Before OLib_initialize CALL\n");
  OLib_initialize();
  printf("After OLib_initialize CALL\n");

  // Load test AES keys into the enclave
  unsigned char inkey[16];
  unsigned char outkey[16];
  read(randfd, inkey, 16);
  read(randfd, outkey, 16);

  // An AES key for constructing and verifying the _plaintext_ of the blocks
  // so that we can ensure that the blocks come out unaltered
  unsigned char datakey[16];
  read(randfd, datakey, 16);
  
  const size_t enc_block_size = 12 + BLOCK_SIZE+ 16;

  unsigned char *inoutbuf = (unsigned char*) malloc(N * enc_block_size);
  unsigned char *outbuf = (unsigned char*) malloc(N * enc_block_size);
  unsigned char *bufend = inoutbuf + N * enc_block_size;
  unsigned char *bufout_end = outbuf + N * enc_block_size;
  size_t buf_len = N * enc_block_size;
 
  printf("Before optimize CALL\n");
  BRN_params params;
  bos_params bosparams; 
  optimize(N, BLOCK_SIZE, -40, &params);


  printf("params.f = %d, params.d = %d, params.B = %d\n", params.f, params.d, params.B);
  //Check contents of params that come out!

  /*
  for (unsigned char *enc_block_ptr = inoutbuf; enc_block_ptr < bufend; enc_block_ptr+=enc_block_size) {
    unsigned char block[BLOCK_SIZE] = {}; // Initializes to zero

    uint64_t rnd;
    read(randfd, (unsigned char*) &rnd, sizeof(rnd));
    //For easier visual debugging:
    rnd = rnd % N;
    memcpy(block, (unsigned char*) &rnd, sizeof(rnd));    
    //printf("%ld, ", rnd);

    // Encrypt the chunk to the enclave      
    read(randfd, enc_block_ptr, 12); // Choose a random IV
    if (BLOCK_SIZE != gcm_encrypt(block, BLOCK_SIZE, NULL, 0, inkey,
            enc_block_ptr, 12, enc_block_ptr+12, enc_block_ptr + 12 + BLOCK_SIZE)) {
        printf("Encryption failed\n");
        break;
    }
  } 
  //printf("\n");
  close(randfd); 
  
  //printf("Before Enclave_loadTestKeys()\n");
  Enclave_loadTestKeys(inkey, outkey);

  bosparams.f = params.f;
  bosparams.n = N;
  bosparams.n_1 = N;
  bosparams.b = pow(params.f, params.d);
  bosparams.Z = ceil(2 * (double)N/(double(bosparams.b)));
  bosparams.d = params.d;
  bosparams.B = params.B;
  bosparams.block_size = BLOCK_SIZE;
  bosparams.evict_start = 0;

  // setup params
  BOS_initialize(&bosparams);
  //double pure_ptime = DecryptAndBOS(inoutbuf, N, enc_block_size, outbuf);
  //printf("pure_ptime = %f\n", pure_ptime);

    #if 1
    // Decrypt the output and do a sanity check
    uint64_t lastkey = 0;
    int numfailed = 0;
    int cnum = 0;
    
    unsigned char* ass_data = NULL;
    unsigned char last_packet_tag[16];
    size_t ass_size = 0;
    for (unsigned char *enc_block_ptr = outbuf; enc_block_ptr < bufout_end; enc_block_ptr += enc_block_size) {
      unsigned char decrypted_block[BLOCK_SIZE];

      #ifdef BOS_OUTSIDE_PRM_STORAGE
        if(enc_block_ptr!=outbuf){
          ass_data = last_packet_tag;
          ass_size = 16;  
        }
      #endif

      ++cnum;
      int gcm_ret = gcm_decrypt(enc_block_ptr+12, BLOCK_SIZE, ass_data, ass_size, enc_block_ptr+12+BLOCK_SIZE, outkey, enc_block_ptr, 12, decrypted_block);
      //printf("outer gcm_ret = %d\n", gcm_ret);
      if (gcm_ret != BLOCK_SIZE) {
          printf("Outer Decryption failed %d/%d\n", ++numfailed, cnum);
      }

      #ifdef BOS_OUTSIDE_PRM_STORAGE    
        memcpy(last_packet_tag, enc_block_ptr+12+chunk_size, 16);
      #endif

      uint64_t blockkey = *(uint64_t *)decrypted_block;
      if (blockkey < lastkey) {
          printf("Unsorted!\n");
      }

      // Note open multi line comment here
      // Check the integrity of the data that came out
      unsigned char should_be_zeroes[sizeof(zeroes)];

      gcm_ret =  gcm_decrypt(chunkbuf+8+12, sizeof(zeroes), chunkbuf, 8,
              chunkbuf+8+12+sizeof(zeroes), datakey, chunkbuf+8, 12, should_be_zeroes);  
      //printf("sizeof(zeroes) = %lu, gcm_ret = %d\n", sizeof(zeroes), gcm_ret);
      bool f1 = sizeof(zeroes) != gcm_ret;
      bool f2 = memcmp(should_be_zeroes, zeroes, sizeof(zeroes));

      if(f1 || f2) {
       //printf("f1 = %d, f2 = %d\n",f1,f2);
      
      //if (sizeof(zeroes) != gcm_decrypt(chunkbuf+8+12, sizeof(zeroes), chunkbuf, 8,
      //        chunkbuf+8+12+sizeof(zeroes), datakey, chunkbuf+8, 12, should_be_zeroes) ||
      //        memcmp(should_be_zeroes, zeroes, sizeof(zeroes))) {
      
          printf("Inner decryption failed\n");
          for(int i=0;i<8;++i) printf("%02x", chunkbuf[i]); printf("\n");
          for(int i=0;i<12;++i) printf("%02x", chunkbuf[8+i]); printf("\n");
          for(int i=0;i<sizeof(zeroes);++i) printf("%02x", chunkbuf[8+12+i]); printf("\n");
          for(int i=0;i<16;++i) printf("%02x", chunkbuf[8+12+sizeof(zeroes)+i]); printf("\n");
      }

      //NOTE close multi line comment here
      lastkey = blockkey;
    }
    #endif
  */
  free(inoutbuf); 
  free(outbuf);
  return 0; 
}
