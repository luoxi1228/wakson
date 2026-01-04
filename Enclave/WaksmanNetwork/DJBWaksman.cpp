#include "DJBWaksman.hpp"

template <typename T> 
void inline reverse_interleaveList(T* list, size_t N, DJB_GCB_state *state){
  //T *temp = new T[N];
  T *temp = (T*)(state->c);
  for(size_t i=0; i<N/2; i++) {
    temp[2*i]= list[i];
  }
  for(size_t i=0; i<N/2; i++) {
    temp[2*i+1]=list[N/2 + i];
  }

  memcpy(list, temp, N *sizeof(T));
}


template <typename T> 
void inline interleaveList_div2(T* list, size_t N, DJB_GCB_state *state){
  //T *temp = new T[N];
  T *temp = (T*)(state->cp);
  for(size_t i=0; i<N/2; i++) {
    temp[i]= list[2*i]/2;
  }
  for(size_t i=0; i<N/2; i++) {
    temp[N/2+i]=list[2*i+1]/2;
  }

  memcpy(list, temp, N *sizeof(T));
}

// Interleave control bits of recursion children.
void inline interleaveControlBits(bool* c_f, bool* c_l, size_t N, size_t gN, DJB_GCB_state *state){
  // skip_value is the number of control bits per layer of the 2*log2(N/2)-1 layers 
  // of one half of the recursion
  int skip_value = N/4;
  int num_bits_per_list = (2 * log2(N/2)-1) * skip_value;
  //printf("num_bits_per_list = %d, skip_value = %d\n", num_bits_per_list, skip_value);
  bool *list = (bool*)(state->cp);
 
  bool *c_start = c_f + gN/2;
  bool *c_ptr = c_start;
  size_t ctr = 0;

  while(ctr < num_bits_per_list){
    list[ctr] = *c_ptr;
    c_ptr++;
    ctr++;
    if(ctr%skip_value==0)
      c_ptr+=(gN/2-skip_value); 
  } 

  c_start = c_start + skip_value;
  c_ptr = c_start;
 
  while(ctr < 2*num_bits_per_list){
    list[ctr] = *c_ptr;
    c_ptr++;
    ctr++;
    if(ctr%skip_value==0)
      c_ptr+=(gN/2-skip_value); 
  } 

  reverse_interleaveList(list, 2 * num_bits_per_list, state);

  c_start = c_f + gN/2;
  c_ptr = c_start;
  ctr=0;

  while(ctr < 2*num_bits_per_list){
    *c_ptr = list[ctr];
    c_ptr++;
    ctr++;
    if(ctr%(2*skip_value)==0)
      c_ptr+=(gN/2- (2*skip_value)); 
  } 
}

//The recursive function to set the control bits of the instantiated Waksman Network
void generateControlBits(uint64_t *pi, uint64_t N, uint64_t gN, bool *controlbits_f, bool *controlbits_l, DJB_GCB_state *state){
  // printf("N = %d\n", N);
  // pi = permutation
  int m = 1;
  while(1<<m < N){
    m+=1;
  }

  // Handle base case:
  // Note that control_bits_f and control_bits_l should point to the same for these.
  if(m==1){
    controlbits_f[0]=pi[0]; 
    return;
  }

  // Initialize p,q, and piinv
  // piinv is initialized to range(N)
  for(size_t i=0; i<N; i++) {
    ((uint64_t*)(state->p))[i]=pi[i^1];
    ((uint64_t*)(state->q))[i]=pi[i]^1;
    (state->piinv)[i]=i;
  }

  // From the paper y = composeinv(a,b) for a,b,y with N elements:
  // returns [y for (x,y) in sorted(zip(b,a))]
  // Translates to the a returned by BitonicSort(b, N, a, NULL, 8, true);
  // but note we get result y as a in place. 
  // So may need a temp buff to store a, if we want to retain former a. 

  // Initialize temp_array_p with pi to compute piinv
  memcpy(state->temp_array_p, pi, N*8);
  // piinv = composeinv(range(n), pi)
  // Sort the temp_array_p (pi) in place to create piinv

  BitonicSort((unsigned char*) state->temp_array_p, N, (unsigned char*) state->piinv, NULL, 8, true); 
  /*
  printf("N=%ld, piinv: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", (state->piinv)[j]);
  }
  printf("]\n");
  */

  // Initialize temp_array_p with p, we'll need this for composeinv(p,q) and composeinv(q,p)
  memcpy(state->temp_array_p, state->p, N*8);
  // Initialize temp_array_q with q, we'll need this for composeinv(p,q) and composeinv(q,p)
  memcpy(state->temp_array_q, state->q, N*8);
  
  //p,q = composeinv(p,q), composeinv(q,p)
  BitonicSort((unsigned char*) state->temp_array_q, N, (unsigned char*) state->p, NULL, 8, true);
  BitonicSort((unsigned char*) state->temp_array_p, N, (unsigned char*) state->q, NULL, 8, true);

  //c = min(x, p[x]) for x in range(n)
  for(size_t i=0; i<N; i++){
    ((uint64_t*)(state->c))[i] = min(i, ((uint64_t*)(state->p))[i]);
  } 
 
  memcpy(state->temp_array_p, state->p, N*8);
  memcpy(state->temp_array_q, state->q, N*8);
  //p,q = composeinv(p,q), composeinv(q,p)
  BitonicSort((unsigned char*) state->temp_array_q, N, (unsigned char*) state->p, NULL, 8, true);
  BitonicSort((unsigned char*) state->temp_array_p, N, (unsigned char*) state->q, NULL, 8, true);

  for(int i=1; i<m-1; i++) {
    //NOTE: We can save the memory for cp, by using temp_array_q as cp, and then
    // refreshing it to q before the next BitonicSort call
    memcpy(state->temp_array_p, state->p, N*8);
    memcpy(state->temp_array_q, state->q, N*8);
    memcpy(state->cp, state->c, N*8);
    // cp, p, q = composeinv(c,q), composeinv(p,q), composeinv(q,p)
    // We do composeinv(c,q) and composeinv(p,q) in one call with this BitonicSort()
    BitonicSort((unsigned char*) state->temp_array_q, N, (unsigned char*) state->cp, (unsigned char*) state->p, 8, true);
    BitonicSort((unsigned char*) state->temp_array_p, N, (unsigned char*) state->q, NULL, 8, true);
    //BitonicSort((unsigned char*) state->p, N, (unsigned char*) state->q, NULL, 8, true);
    for(size_t j=0; j<N; j++){
      ((uint64_t*)(state->c))[j]= min(((uint64_t*)(state->c))[j],((uint64_t*)(state->cp))[j]);
    }  
  }

  /*
  printf("N=%ld, c: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", ((uint64_t*)(state->c))[j]);
  }
  printf("]\n");
  */

  // Computing f and F in one loop:
  for(int j=0; j<N/2; j++){
    //Convert to shift operation
    //state->p is f
    ((bool*)(state->p))[j] = ((uint64_t*)(state->c))[2*j]%2;

    (state->F)[2*j] = (2*j) ^ (((bool*)(state->p))[j]);
    (state->F)[2*j+1] = (2*j+1) ^ (((bool*)(state->p))[j]);

  }

  /*
  printf("N=%ld, piinv: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", (state->piinv)[j]);
  }
  printf("]\n");

  printf("N=%ld, F: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", (state->F)[j]);
  }
  printf("]\n");
  */

  // Fpi = composeinv(F, piinv)
  BitonicSort((unsigned char*) state->piinv, N, (unsigned char*) state->F, NULL, 8, true);
  // F is now Fpi  

  /*
  printf("N=%ld, Fpi: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", (state->F)[j]);
  }
  printf("]\n");
  */   

  for(int j=0; j<N/2; j++){
    //Convert to shift operation
    // state->q is l
    ((bool*)(state->q))[j] = (state->F)[2*j]%2;
    
    ((size_t*)(state->c))[2*j] = (2*j) ^ (((bool*)(state->q))[j]);
    ((size_t*)(state->c))[2*j+1] = (2*j+1) ^ (((bool*)(state->q))[j]);
  }

  /*
  printf("N=%ld, L: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", ((size_t*)(state->c))[j]);
  }
  printf("]\n");
  */

  // M = composeinv(Fpi, L)
  BitonicSort((unsigned char*) state->c, N, (unsigned char*) state->F, NULL, 8, true);
  // F is now M, since Fpi was F

  /*
  printf("N=%ld, M: [", N);
  for(int j = 0; j < N; j++) {
    printf("%d, ", (state->F)[j]);
  }
  printf("]\n");
  */

  interleaveList_div2(state->F, N, state);
  // F = Fpi = M is now interleaved such that all the even indexes are grouped to
  // the first half of the array, and the odd ones in the second half.
 
  // Set the control bits of f and l  
  for(int j=0; j<N/2; j++){
    controlbits_f[j] = ((bool*)(state->p))[j];
    controlbits_l[j] = ((bool*)(state->q))[j];
  }

  size_t *pi_next = state->F; 
  state->F = state->F + gN;
  generateControlBits(pi_next, N/2, gN, controlbits_f + (gN/2), controlbits_l - gN/2, state);
  state->F = state->F + N/2;
  generateControlBits(pi_next + N/2, N/2, gN, controlbits_f + gN/2 + N/4, controlbits_l -gN/2 + N/4, state);
  state->F = state->F - (gN + (N/2));

  // Interleave the returned control bits from the recursion!
  if(m>1){
    interleaveControlBits(controlbits_f, controlbits_l, N, gN, state);
  }  
}


void generateControlBits(uint64_t *pi, uint64_t N, uint64_t gN, bool *controlbits_f, bool *controlbits_l){
  DJB_GCB_state *state = new DJB_GCB_state(N);
  generateControlBits(pi, N, gN, controlbits_f, controlbits_l, state);
  delete state;
}

void DJBWaksmanShuffle(unsigned char *buffer, size_t N, size_t block_size, enc_ret *ret) { 
  // Note permutations use 8 byte values for indices in DJB Waksman modes
  // and 4 bytes in OA Waksman modes. 
  uint64_t *random_permutation;  
  try {
    random_permutation = new uint64_t[N];
  } catch (std::bad_alloc&){
    printf("Allocating memory failed in DJBWaksmanShuffle\n");
  }

  long t1, t2;
  ocall_clock(&t1);
  generateRandomPermutation(N, random_permutation);
  ocall_clock(&t2);
  double time_rp_ms = ((double)(t2-t1))/1000.0;
  ret->gen_perm_time = time_rp_ms;
  ret->OSWAP_gp=OSWAP_COUNTER;
  OSWAP_COUNTER=0; 

  #ifdef TEST_WN_DJB
    uint64_t *correct_permuted_keys = new uint64_t[N];
    for(size_t i=0; i<N; i++) {
      uint64_t buffer_key = *((uint64_t*)(buffer + (block_size * random_permutation[i])));
      correct_permuted_keys[i] = buffer_key;
      //printf("correct_permuted_keys = %ld, buffer_key = %ld\n",correct_permuted_keys[i], buffer_key);
    }
  #endif  

  int logn = calculatelog2(N); 
  size_t num_switches = (N/2) * (2*logn-1); 
  bool *controlbits = new bool[num_switches];
  ocall_clock(&t1);
  generateControlBits(random_permutation, N, N, controlbits, controlbits + (N*(logn-1)));
  ocall_clock(&t2);
  double time_gcb_ms = ((double)(t2-t1))/1000.0;
  ret->control_bits_time = time_gcb_ms;
  ret->OSWAP_cb=OSWAP_COUNTER;
  OSWAP_COUNTER=0;

  // Testing generated controlBits:
  /*
  printf("\nGenerated Control Bits :\n");
  for(size_t i=0; i<num_switches; i++){
    if(i%(N/2)==0)
      printf("\n");
    printf("%d, ", controlbits[i]);
  }
  printf("\n");
  */

  // The old non-recursive, layer-by-layer Apply Permutation is commented out below. 
  bool *controlbits_b = controlbits+((2*logn-2)*N/2);
  ocall_clock(&t1);
  if(block_size==8){
    applyPermutation<OSWAP_8>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
    //applyPermutation<OSWAP_8>(controlbits, num_switches, buffer, N, block_size);
  } else if(block_size%16==0){
    applyPermutation<OSWAP_16X>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
    //applyPermutation<OSWAP_16X>(controlbits, num_switches, buffer, N, block_size);
  } else {
    applyPermutation<OSWAP_8_16X>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
    //applyPermutation<OSWAP_8_16X>(controlbits, num_switches, buffer, N, block_size);
  }
  ocall_clock(&t2);
  double time_ap_ms = ((double)(t2-t1))/1000.0;
  ret->apply_perm_time = time_ap_ms;
  ret->OSWAP_ap=OSWAP_COUNTER;
  
  #ifdef TEST_WN_DJB
    unsigned char *buffer_ptr = buffer;
    bool failed = false;
    for(size_t i=0; i<N; i++) {
      uint64_t buffer_key = *((uint64_t*)(buffer_ptr));
      //printf("correct_permuted_keys = %ld, buffer_key = %ld\n",correct_permuted_keys[i], buffer_key);
      if(correct_permuted_keys[i]!=buffer_key) {
        printf("TEST_WN_DJB: Shuffle Correctness Failed\n");
        failed = true;
        break;
      }
      buffer_ptr+=block_size;
    }
    //if(!failed)
      //printf("TEST_WN_DJB: Shuffle Correctness SUCCESS!\n");
    delete []correct_permuted_keys;
  #endif  

  delete []random_permutation;
  delete []controlbits;
}

void DecryptAndDJBWaksmanShuffle(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret) {
  long t1, t2;

  // Decrypt buffer to decrypted_buffer
  unsigned char *decrypted_buffer = NULL;
  size_t decrypted_block_size = decryptBuffer(encrypted_buffer, N, encrypted_block_size,
    &decrypted_buffer);
  ocall_clock(&t1);

  PRB_pool_init(1);
  DJBWaksmanShuffle(decrypted_buffer, N, decrypted_block_size, ret);
  ocall_clock(&t2);

  // Encrypt buffer to result_buffer
  encryptBuffer(decrypted_buffer, N, decrypted_block_size, result_buffer);
  PRB_pool_shutdown();

  double compare_ms = ((double)(t2-t1))/1000.0;
  //double encryption_ms = ((double)(t3-t2))/1000.0;

  free(decrypted_buffer);
  ret->ptime = compare_ms;
  ret->OSWAP_count = OSWAP_COUNTER;
  return;
}

void DJBWaksmanSort(unsigned char *buffer, size_t N, size_t block_size, enc_ret *ret) {
  long t1, t2;
  size_t *permutation; 
  try {
    permutation = new size_t [N];
  } catch (std::bad_alloc&) {
    printf("Allocating memory failed in DJBWaksmanSort\n");
  }
  unsigned char *buffer_ptr = buffer;

  ocall_clock(&t1);
  generateSortPermutation_DJB(N, buffer, block_size, permutation);
  ocall_clock(&t2);
  double gen_perm_ms = ((double)(t2-t1))/1000.0;
  ret->gen_perm_time = gen_perm_ms;
  ret->OSWAP_gp=OSWAP_COUNTER;
  OSWAP_COUNTER=0;
 
  int logn = calculatelog2(N);
  size_t num_switches = (N/2) * (2*logn-1);
  bool *controlbits = new bool[num_switches];
  ocall_clock(&t1);
  generateControlBits(permutation, N, N, controlbits, controlbits + (N*(logn-1)));
  ocall_clock(&t2);
  double time_gcb_ms = ((double)(t2-t1))/1000.0;
  ret->control_bits_time = time_gcb_ms;
  ret->OSWAP_cb=OSWAP_COUNTER;
  OSWAP_COUNTER=0;

  // Testing generated controlBits: 
  /*
  printf("\nGenerated Control Bits :\n");
  for(size_t i=0; i<num_switches; i++){
    if(i%(N/2)==0)
      printf("\n");
    printf("%d, ", controlbits[i]);
  }
  printf("\n");
  */

  ocall_clock(&t1);
  bool *controlbits_b = controlbits+((2*logn-2)*N/2);
 
  if(block_size==8){
    applyPermutation<OSWAP_8>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
  } else if(block_size%16==0){
    applyPermutation<OSWAP_16X>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
  } else {
    applyPermutation<OSWAP_8_16X>(controlbits, controlbits_b, buffer, N, N, block_size, 0);
  }
  ret->OSWAP_ap=OSWAP_COUNTER;

  ocall_clock(&t2);
  double time_ap_ms = ((double)(t2-t1))/1000.0;
  ret->apply_perm_time = time_ap_ms;
  
  delete []controlbits;
  delete []permutation;
}

void DecryptAndDJBWaksmanSort(unsigned char *encrypted_buffer, size_t N, size_t encrypted_block_size, unsigned char *result_buffer, enc_ret *ret) {
  // Decrypt buffer to decrypted_buffer
  long t1, t2;
  unsigned char *decrypted_buffer = NULL;
  size_t decrypted_block_size = decryptBuffer(encrypted_buffer, N, encrypted_block_size,
    &decrypted_buffer);
  ocall_clock(&t1);

  PRB_pool_init(1);
  DJBWaksmanSort(decrypted_buffer, N, decrypted_block_size, ret);
  ocall_clock(&t2);

  // Encrypt buffer to result_buffer
  encryptBuffer(decrypted_buffer, N, decrypted_block_size, result_buffer);
  PRB_pool_shutdown();

  double compare_ms = ((double)(t2-t1))/1000.0;
  //double encryption_ms = ((double)(t3-t2))/1000.0;

  free(decrypted_buffer);
  ret->ptime = compare_ms;
  ret->OSWAP_count = OSWAP_COUNTER;
  return;

}
