#!/usr/bin/python3

import subprocess
import os
import math
import sys

###############################################################################

# CONFIGS to SET:

RESULTS_FOLDER = "../RESULTS"

MODE = [6] #[10, 12, 13, 14, 15, 17, 2, 3, 5, 7, 9]

N = [65536, 262144, 1048576]
BLOCK_SIZE = [64, 256, 1024]

REPEAT = 5

# Optimization target for BORP/BOS:
# 1 = Optimize Total time
# 2 = Optimize Phase 2 time
OPT_TARGET = 2

###############################################################################

# List of Modes:
'''
    Shuffles:

    1 = RecursiveShuffle_M1 (with Goodrich's OP_TightCompact)
    2 = ORShuffle
    3 = BORPStream (Bucket Oblivious Random Permutation Stream) (lambda = -80)
    4 = OddEvenMergeSort Shuffle
    5 = BitonicSort Shuffle
    6 = ButterflyNetwork Compaction
    7 = WaksShuffle
    8 = Bucket Oblivious Random Permutation (with TC)
    9 = Nassimi-Sahni Waksman Shuffle

    ---------------------------------------------------------------------------
    Sorts:

    10 = Nassimi-Sahni Waksman Sort
    11 = Oblivious Sort (Sorting Network - OddEvenMergeSort)
    12 = Oblivious Sort (Sorting Network - BitonicSort)
    13 = Oblivious Sort (BORPStream + Quicksort, lambda = -80)
    14 = Oblivious Sort (ORShuffle + Quicksort)
    15 = WaksShuffle (Mode 7) + Quicksort
    17 = WaksSort
'''

###############################################################################

# Some useful common parameters for quickly setting up experiments:

'''
#BLOCK_SIZE = [8, 16, 32, 64, 128, 256, 304, 512, 768, 1024]
#N = [1000, 5000, 10000, 20000, 50000, 60000, 70000, 80000, 90000, 100000]
#N = [100000, 250000, 500000, 750000, 1000000, 1250000, 1500000, 1750000, 2000000, 
#     3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000]
'''

'''
#Powers of 2 N:
N = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
     1048576, 2097152, 4194304, 8388608, 16777216]
'''

'''
#Powers of 2 N (with 2 points in between x1.26, x1.59)
N = [1024, 1290, 1638,  2048, 2580, 3256, 4096, 5161, 6513, 8192, 10322, 13025,
     16384, 20644, 26051, 32768, 41288, 52101, 65536, 82575, 104202, 131072, 
     165151, 208404, 262144, 330301, 416808, 524288, 660603, 833618, 1048576,
     1321206, 1667236, 2097152, 2642412, 3334472, 4194304, 5284823, 6668943, 
     8388608, 10569646, 13337887, 16777216]
'''

#N = [32768, 131072, 1048576]
#BLOCK_SIZE = [64, 128, 256, 512, 1024]

###############################################################################

def calculateZ_prime(N, f, d, B):
  denom = math.pow(f,d)
  Z_prime = math.ceil(float(2*N)/float(denom)) + d*B
  #print ("Z_prime = %d, for f = %d, d = %d, B = %d" %(Z_prime, f, d, B))
  return Z_prime
if(os.path.exists(RESULTS_FOLDER)==False):
  os.mkdir(RESULTS_FOLDER)

# 向上取最近的2的幂
def pow2_gt(x):  
  log2 = math.log(x, 2.0)
  log2_ceil = math.ceil(log2)
  return (math.pow(2, log2_ceil))

BASE_HEAP = 1000000
# Enclave.config.xml generator:
# We need to generate a config file with heap size: N*block_size + N + base
#   - N * block_size : required to hold decrypted buffer in memery in enclave
#   - N : for TC (Goodrich's) which needs the offset array
#   - base : for anyother small instance of heap usage. Setting this to 10KB for now
#
# This value is then rounded up to the next multiple of 4096 (4KiB) for matching page requirements of SGX

def generateConfigFile(MODE, N, block_size):
  heap_memory = (N * block_size) + (2 * N * 8) + BASE_HEAP
  
  #Accounting for random_bytes_buffer memory usage, now that we grab random_bytes in bulk upfront.
  if(MODE==1 or MODE==2 or MODE==12 or MODE==14):
    #MAKE sure this additional constant heap memory matches the RS_RB_BUFFER_SIZE from CONFIG 
    heap_memory+=1000000

  if (MODE == 6):
    logN = math.log(N, 2)
    # 1) 控制位/交换位：大约 N/2 * logN 个开关，按 1 byte/4 byte 计都行，先保守
    heap_memory += int(N * logN * 32)   # 经验上先给 32NlogN bytes
    # 2) 如果 butterfly 内部要额外 buffer（常见），给 1-2 倍明文 buffer 的空间
    heap_memory += int(2 * N * block_size)
    heap_memory += 1000000


  if(MODE in {7,15,17}):
    logN = math.log(N,2)
    #heap_memory+=(N*4) #For the random_permutation array
    #For the switches (4 per wire, 1 for in & out switch (n switches in total) = 5, 16 from pointers to subnw)
    # Updating it:
    heap_memory+=(N*logN * 21)    
    #heap_memory+=logN * (4 + N * 21)   
    #For the setControlBits (from forward_perm, unselected_cnt, input_sw)
    heap_memory+=(2*N*48) #For 2x forward and rev_perms
    heap_memory+=(N*12) #For unselected_cnt, input_sw
    heap_memory+=1000000

  if(MODE==8):
    N_prime = pow2_gt(N)
    print(N_prime)
    # For all the buckets
    heap_memory+=(N_prime * 2 * (block_size+8))
    heap_memory+=(2 * 256 * (block_size+8))
    # For TC:
    heap_memory+=(2 * N * 8)
    #MAKE sure this additional constant heap memory matches the RS_RB_BUFFER_SIZE from CONFIG 
    heap_memory+=1000000

  if(MODE in {9,10}):
    #For controlbits: N/2 * 2*logn-1
    logN = math.log(N,2)
    heap_memory+= (N * logN)
    # generateControlBits(): for the additional lists required:
    # p+q+tempp+tempq+piinv = 5 * 8 * N, cp, c = 2 * N * logN, F = 8 * N * logN
    if(logN<8):
      heap_memory+=((5 * 8 * N) + (2 *  N * 8) + (8 * N * logN))
    else:
      heap_memory+=((5 * 8 * N) + (2 *  N * logN) + (8 * N * logN))
    # For the RS:
    heap_memory+=(9*N)
    #MAKE sure this additional constant heap memory matches the RS_RB_BUFFER_SIZE from CONFIG 
    heap_memory+=1000000  

  if(MODE==61):
    logN = math.log(N, 2)
    heap_memory+=(8 * N * logN)
    heap_memory+=N

  if (MODE in {3,13}):
   
    # Updating <f, d, s> selection with new failure probability calculator.
    # Currently just doing it for lmbda = -80.
    
    # Old incorrect Markov calculator:
    #f_opt, d_opt, B_opt = returnBOSParams(N, block_size, lmbda, OPT_TARGET)
    f_opt = 2
    d_opt = 1
    B_opt = 29

    if(OPT_TARGET==1):
      if(block_size == 8):
        if(N <= 4096):
          B_opt = 29
        elif(N <= 16384):
          B_opt = 30
        elif(N <= 262144):
          B_opt = 31
        elif(N <= 2097152):
          B_opt = 32
        elif(N <= 8388608):
          B_opt = 33
        else:
          f_opt = 4
          B_opt = 48
      else:
        #NOTE: We implicitly assume if block_size !=8, N = 1048576.
        if(block_size <= 32):
          B_opt = 32
        elif(block_size <= 128):
          f_opt = 4
          B_opt = 46
        elif(block_size <= 512):
          f_opt = 4
          d_opt = 2
          B_opt = 47
        else:
          f_opt = 4
          d_opt = 3
          B_opt = 47

    elif(OPT_TARGET==2):
      if(block_size==8):
        if(N <= 2048):
          B_opt = 29
        if(N <= 4096):
          d_opt = 2
          B_opt = 30
        if(N <= 8192):
          f_opt = 4
          B_opt = 43
        if(N < 16384):
          f_opt = 4
          d_opt = 2
          B_opt = 44
        if(N <= 65536):
          f_opt = 4
          d_opt = 2
          B_opt = 45
        if(N <= 262144):
          f_opt = 4
          d_opt = 3
          B_opt = 46
        if(N <= 524288):
          f_opt = 4
          d_opt = 3
          B_opt = 47
        if(N <= 1048576):
          f_opt = 4
          d_opt = 4
          B_opt = 47
        if(N <= 4194304):
          f_opt = 4
          d_opt = 4
          B_opt = 48
        else:
          f_opt = 4
          d_opt = 5
          B_opt = 49
      else: 
        #NOTE: We implicitly assume if block_size !=8, N = 1048576.
        f_opt = 4
        d_opt = 4
        B_opt = 47

    b_opt = math.pow(f_opt, d_opt-1)
    print("f_opt = %d, d_opt = %d, B_opt = %d"%(f_opt, d_opt, B_opt))

    #For the decrypted buffer 
    heap_memory += (N*block_size)
    #For all the MSN and FN buffers:
    # + 100 to account for the other internal variables of all the MSN/FNs
    heap_memory += ((B_opt * (block_size +16) + 100) *d_opt * b_opt)
    #For the outbuf with all the dummies + additional packets from flush buffers:
    Zprime = calculateZ_prime(N, f_opt, d_opt, B_opt)
    heap_memory+= (Zprime * block_size * (f_opt**d_opt))  
    #For the added RS at the end:
    heap_memory+=(Zprime*math.log(Zprime,2)*8)
    # 1M for RS's local PRB, 1M for PRB for labels
    heap_memory+=2000000
    #print("heap_memory = %d"%(heap_memory))
    #print(f_opt, d_opt, B_opt, b_opt)



  # 向上对齐到 4KiB 页的倍数
  num_pages = int(heap_memory/4096)
  if(int(heap_memory)%4096!=0):
    heap_memory = (num_pages+1) * 4096
  heap_memory_hex = hex(int(heap_memory))
  heap_line = "  <HeapMaxSize>"+str(heap_memory_hex)+"</HeapMaxSize>\n"
  #print(heap_line)
  config_file = open("../Enclave/Enclave.config.xml","r")
  lines = config_file.readlines()
  config_file.close()
  lines[4] = heap_line
  config_file = open("../Enclave/Enclave.config.xml","w")
  config_file.writelines(lines)
  config_file.close()

  #print("heap_memory = %d"%(heap_memory))
  #print("heap_memory = " + hex(heap_memory))
  #print(heap_line)
  #skip_exp = heap_memory > SKIP_LIMIT
  skip_exp = False
 
  if(MODE==3 or MODE==13):
    return skip_exp, f_opt, d_opt, B_opt
  else:
    return skip_exp, 0, 0, 0


# Configure Experiment Parameters:
for m in MODE:
  skipped = []
  for b in BLOCK_SIZE:
    ttime = []
    ptime = []      
    lltime = []
    ttime_stddev = []
    ptime_stddev = []
    lltime_stddev = []
    oswap_count = []
    set_cb = []
    set_cb_stddev = []
    apply_perm = []
    apply_perm_stddev = [] 
    oswap_count_gp = []
    oswap_count_cb = []
    oswap_count_ap = []
    qtime = []
    qtime_stddev = []

    for n in N:
      print("Running experiment m = %d, b = %d, n = %d" % (m,b,n))
      skip, f, d, B = generateConfigFile(m, n, b)
      if(skip):
        print("Skipping experiment (heap_memory request > SKIP_LIMIT)")
        exp = {n,b}
        skipped.append(exp)
        continue

      p = subprocess.run(["make", "-C", "../"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

      # Run script_application with the corresponding m, n, and b, and log output
      # Parse returned output
      if(m in {3,13}):
        p = subprocess.run(["./script_application", str(m), str(n), str(b), str(REPEAT), str(f), str(d), str(B)], stdout=subprocess.PIPE)  
      else:
        p = subprocess.run(["./script_application", str(m), str(n), str(b), str(REPEAT)], stdout=subprocess.PIPE)
        # p = subprocess.run(["./script_application", str(m), str(n), str(b), str(REPEAT)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

      # if p.returncode != 0:
      #     print(f"CRITICAL ERROR: Application failed with return code {p.returncode}")
      #     print(f"STDERR content:\n{p.stderr.decode('utf-8')}")

      out_lines = (p.stdout.decode("utf-8").split('\n'))
      print(out_lines)
      if(m in {7,9,10,17}):
        out_lines = out_lines[-15:-1]
      elif(m==15):
        out_lines = out_lines[-17:-1]
      elif(m==13):
        out_lines = out_lines[-12:-1]
      elif(m==3 or m==14):
        out_lines = out_lines[-10:-1]
      elif(m == 6):
        out_lines = out_lines[-14:-1]
      else:
        out_lines = out_lines[-8:-1]
      print("Out_lines:")
      print(out_lines)
      '''
      flag_incorrect_num_lines = False
      if(not(len(out_lines)==7 or len(out_lines)==13)):
        print(len(out_lines))
        flag_incorrect_num_lines = True
      '''

      flag_notFloat = False
      for line in out_lines:
        try:
          float(line.strip())
        except ValueError:
          print("Line with value error is:")
          print(line)
          flag_notFloat = True
          break  
    
      #if(flag_incorrect_num_lines or flag_notFloat):
      if(flag_notFloat):
        print("Receieved unexpected output, this experiment run has failed. ONE MUST DEBUG!")
        continue
      else:
        ttime.append(out_lines[0].strip())
        ttime_stddev.append(out_lines[1].strip()) 
        ptime.append(out_lines[2].strip())
        ptime_stddev.append(out_lines[3].strip())
        oswap_count.append(int(out_lines[4].strip()))
        lltime.append(out_lines[5].strip())
        lltime_stddev.append(out_lines[6].strip())
        
        if(m==3):
          oswap_count_cb.append(int(out_lines[7].strip()))
          oswap_count_ap.append(int(out_lines[8].strip()))

        if(m==14):
          qtime.append(out_lines[7].strip())      
          qtime_stddev.append(out_lines[8].strip())      
        
        if(m==13):
          oswap_count_cb.append(int(out_lines[7].strip()))
          oswap_count_ap.append(int(out_lines[8].strip()))
          qtime.append(out_lines[9].strip())      
          qtime_stddev.append(out_lines[10].strip())      
  
        if(m in {7,9,10,15,17}):
          set_cb.append(out_lines[7].strip())
          set_cb_stddev.append(out_lines[8].strip()) 
          apply_perm.append(out_lines[9].strip())
          apply_perm_stddev.append(out_lines[10].strip()) 
          oswap_count_gp.append(int(out_lines[11].strip()))
          oswap_count_cb.append(int(out_lines[12].strip()))
          oswap_count_ap.append(int(out_lines[13].strip()))
          if(m==15):
            qtime.append(out_lines[14].strip())      
            qtime_stddev.append(out_lines[15].strip())      

        if(m==6):
          set_cb.append(out_lines[7])
          set_cb_stddev.append(out_lines[8])
          apply_perm.append(out_lines[9])
          apply_perm_stddev.append(out_lines[10])
          oswap_count_cb.append(int(out_lines[11]))
          oswap_count_ap.append(int(out_lines[12]))


    # Log ttime, ptime and their stddevs into log_file_n
    log_file_n = RESULTS_FOLDER+"/"+str(m)+"_"+str(b)+'.lg'         
    log_file = open(log_file_n, "a")
    for i in range(len(N)):
      if({N[i],b} not in skipped):
        if(m in {7,9,10,17}):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+","+str(set_cb[i]) +","+str(set_cb_stddev[i])+","+str(apply_perm[i]) +","+str(apply_perm_stddev[i])+","+str(oswap_count_gp[i])+","+str(oswap_count_cb[i])+","+str(oswap_count_ap[i])+"\n"
        elif(m==15):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+","+str(set_cb[i]) +","+str(set_cb_stddev[i])+","+str(apply_perm[i]) +","+str(apply_perm_stddev[i])+","+str(oswap_count_gp[i])+","+str(oswap_count_cb[i])+","+str(oswap_count_ap[i])+","+str(qtime[i])+","+str(qtime_stddev[i])+"\n"
        elif(m==13):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i]) +","+str(oswap_count_cb[i])+","+str(oswap_count_ap[i])+","+str(qtime[i])+","+str(qtime_stddev[i])+"\n"
        elif(m==14):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+","+str(qtime[i])+","+str(qtime_stddev[i])+"\n"
        elif(m==3):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+","+str(oswap_count_cb[i])+","+str(oswap_count_ap[i])+"\n"
        elif(m==6):
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+","+str(set_cb[i]) +","+str(set_cb_stddev[i])+","+str(apply_perm[i]) +","+str(apply_perm_stddev[i])+","+str(oswap_count_cb[i])+","+str(oswap_count_ap[i])+"\n"
        else:
          log_line = str(N[i])+","+str(ttime[i])+","+str(ttime_stddev[i])+","+str(ptime[i])+","+str(ptime_stddev[i])+","+str(lltime[i])+","+str(lltime_stddev[i])+","+str(oswap_count[i])+"\n"
      else:
        continue
      log_file.write(log_line)
    log_file.close()
