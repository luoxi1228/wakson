#!/usr/bin/python3

import sys
import os
import subprocess
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import math

###############################################################################

# CONFIGS:

# PLOT_MODES:
# 1: Shuffle Algorithms, x-axis: N, y-axis : time (in s) 
# 2: Shuffle Algorithms, x-axis: block_size, y-axis : time (in s) 
# 3: Shuffle Algorithms, x-axis: N, y-axis : OSWAP count

# -----------------------------------------------------------------------------

# 4: Sort Algorithms, x-axis: N, y-axis: time (in s)
# 5: Sort Algorithms, x-axis: block_size, y-axis: time (in s)
# 6: Sort Algorithms, x-axis: N, y-axis: OSWAP count

#Set N and BLOCK_SIZES correctly under the next section of configs 
PLOT_MODE = 4

RESULTS_FOLDER = "Our_Data"

# Set scale of axis:
X_AXIS_LOGSCALE = True
Y_AXIS_LOGSCALE = True

# For modes that use QSort at the end
PLOT_QSORT_TIME = False

# For Waksman modes:
PLOT_GEN_PERM = False
PLOT_CONTROL_BITS = False
PLOT_APPLY_PERM = False
PLOT_OFFLINE_TIME = False

# For Waksman and BORP modes:
PLOT_ONLINE_TIME = True
# For BORP modes ONLINE_TIME is Phase 2 time

PLOT_TIME_IN_MS = False

###############################################################################

if(len(sys.argv)>1):
  RESULTS_FOLDER = (sys.argv[1])
  PLOT_MODE = int(sys.argv[2])

###############################################################################


# Select other parameters depending on CONFIGS

# 1) Assumes PLOT_MODEs 1 and 3 will only have a BLOCK_SIZE array of len 1.
# 2) Assumes PLOT_MODEs 2 and 4 will only have a N array of len 1.
# i.e. PLOT_MODEs 1 and 3 are meant to varying N eon the x-axis against time 
# taken on y-axis and similary 2 and 4 are meant to have varying block_size 
# on x-axis againts varying time on y-axis

flag_BV1 = True
if(not(os.path.isfile(RESULTS_FOLDER+"/"+"33_64.lg"))):
  flag_BV1 = False

if (PLOT_MODE == 1):
  MODE = [9, 3, 33, 2, 7]
  N = [32768, 65536, 131072, 262144, 524288, 1048576]
  #N = [32768, 65536, 131072, 262144]
  # There should only be 1 value in BLOCK_SIZE for this PLOT_MODE
  BLOCK_SIZE = [4096]
  LABEL_STRINGS = ["Nassimi-Sahni", "BORPStream (V2)", "BORPStream (V1)", "ORShuffle", "WaksShuffle" ]
  bs_val = str(BLOCK_SIZE[0])
  save_location = RESULTS_FOLDER+ '/Shuffles_PM1_'+bs_val+'.png'

elif (PLOT_MODE == 2):
  MODE = [9, 3, 33, 2, 7]
  # There should only be 1 value in N for this PLOT_MODE
  N = [1048576]
  BLOCK_SIZE = [64, 128, 256, 512, 1024,
      1152, 1280, 1408, 1536, 1664, 1792, 1920,
      2048, 3072, 4096]
  LABEL_STRINGS = ["Nassimi-Sahni", "BORPStream (V2)" , "BORPStream (V1)", "ORShuffle", "WaksShuffle"]
  if(not(flag_BV1)):
    MODE = [9, 3, 2, 7]
    LABEL_STRINGS = ["Nassimi-Sahni", "BORPStream (V2)", "ORShuffle", "WaksShuffle"]
  LABELS_Y = [64, 128, 256, 512, 1024, 2048, 4096]
  n_val = str(N[0])
  save_location = RESULTS_FOLDER+ '/Shuffles_PM2_'+n_val+'.png'

elif (PLOT_MODE == 3):
  # For OSWAP counts
  print("Not setup yet")
  exit()

elif (PLOT_MODE == 4):
  MODE = [10, 43, 17, 14, 12, 13, 15]
  LABEL_STRINGS = ["Nassimi-Sahni", "BORPStream (V1) + QS", "WaksSort", "ORShuffle + QS", "Bitonic Sort", "BORPStream (V2)\n + QS", "WaksShuffle + QS"] 
  save_location = RESULTS_FOLDER+ '/Sorts_PM4_1024.png'
  #N = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]
  N = [32768, 65536, 131072, 262144, 524288, 1048576]
  # There should only be 1 value in BLOCK_SIZE for this PLOT_MODE
  BLOCK_SIZE = [4096]
  bs_val = str(BLOCK_SIZE[0])
  save_location = RESULTS_FOLDER+ '/Sorts_PM4_'+bs_val+'.png'

elif (PLOT_MODE == 5):
  MODE = [10, 43, 17, 14, 12, 13, 15]
  LABEL_STRINGS = ["Nassimi-Sahni", "BORPStream (V1) + QS", "WaksSort", "ORShuffle + QS", "Bitonic Sort", "BORPStream (V2)\n + QS", "WaksShuffle + QS"]

  if(not(flag_BV1)):
    MODE = [10, 17, 14, 12, 13, 15]
    LABEL_STRINGS = ["Nassimi-Sahni", "WaksSort", "ORShuffle + QS", "Bitonic Sort", "BORPStream (V2)\n + QS", "WaksShuffle + QS"]
  # There should only be 1 value in N for this PLOT_MODE
  N = [32768]
  BLOCK_SIZE = [64, 128, 256, 512, 1024,
      1152, 1280, 1408, 1536, 1664, 1792, 1920,
      2048, 3072, 4096]
  #LABELS_Y = [64, 128, 256, 512, 1024, 2048, 3072, 4096]
  LABELS_Y = [64, 128, 256, 512, 1024, 2048, 4096]
  n_val = str(N[0])
  #save_location = RESULTS_FOLDER+ '/Sorts_PM5_'+n_val+'_total_linlin.png'
  save_location = RESULTS_FOLDER+ '/Sorts_PM5_'+n_val+'_total_lglg.png'

elif (PLOT_MODE == 6):
  # For OSWAP counts
  MODE = [10, 12, 13, 14, 15, 17]
  LABEL_STRINGS = ["DJB Sort", "Bitonic Sort", "BORPStream + QS", "ORShuffle + QS", "OA (Shuffle) + QS", "WaksOn + WaksOff"]
  save_location = RESULTS_FOLDER+ '/Sorts_PM6.png'

  N = [32768, 65536, 131072, 262144, 524288, 1048576]
  BLOCK_SIZE = [1024]

time_divisor = 1000.0
if(PLOT_TIME_IN_MS):
  time_divisor = 1.0
  
###############################################################################

POINT_STYLES = [".","v","s","D","^","d",">","x","*","o","D"]
LINE_STYLES = ['solid', 'dotted', 'dashed', 'dashdot']

PS_SIZE = len(POINT_STYLES)
LS_SIZE = len(LINE_STYLES)
# Colors range from C0 to C9
C_SIZE = 10
#COLOR_MAP_TOTAL 
CMT = {2:1, 3:4, 7:5, 9:6, 10:6, 5:8, 12:8, 13:4, 14:1, 15:5, 17:9, 33:0, 43:0}
#COLOR_MAP_ONLINE 
CMO = {3:3, 7:2, 9:7, 15:2, 13:3}


###############################################################################

# List of Modes:
'''
    Shuffles:

    1 = RecursiveShuffle_M1 (with Goodrich's OP_TightCompact)
    2 = ORShuffle
    3 = BORPStream (Bucket Oblivious Random Permutation Stream) (lambda = -80)
    4 = OddEvenMergeSort Shuffle
    5 = BitonicSort Shuffle
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

# Plot font size manipulators:
title_font_size = 22
legend_font_size = 18
tick_font_size = 18
axis_label_font_size = 22

###############################################################################

TIME_PLOT_MODES = {1, 2, 4, 5}
OSWAP_PLOT_MODES = {3, 6}

ONLY_PTIME_MODES = [1, 2, 4, 5, 33, 43]
BORP_MODES = [3, 13]
WAKSMAN_MODES = [7, 9, 10, 15, 17]
WAKSMAN_SHUFFLE_MODES = [7, 9, 15]
WAKSMAN_SORT_MODE = [10, 17]
QSORT_MODES = [13, 14, 15]

###############################################################################

# Functions:

def convertTo2PowLabels(list_2pows):
  labels_list = []
  for i in list_2pows:
    label = "$2^{"
    log_i = float("{:.1f}".format(math.log2(i)))
    if(log_i.is_integer()):
      log_i = int(math.log2(i))
    label=label+str(log_i)+"}$"
    labels_list.append(label)
  return labels_list
#ax.set_xticklabels([r'$2^6$', r'$2^7$', r'$2^8$', r'$2^9$', r'$2^{10}$'])

###############################################################################


label_ctr=0
color_ctr=0
ps_ctr=0
ls_ctr=0

plt.figure(figsize=(8,5))
# Use the dimensions below for the wider Sort_PM4_4096 graph
#plt.figure(figsize=(9,5))
plt.grid(True)

for m in MODE: 
  ttime = []
  ptime = []
  ttime_stddev= []
  ptime_stddev= []
  OSWAP_count = []
  lltime = []
  lltime_stddev = []
  time_gp = []
  time_cb = []
  time_ap = []
  OSWAP_count_gp = []
  OSWAP_count_cb = []
  OSWAP_count_ap = []
  qsort_time = []
  qsort_time_stddev = []
  online_time = []
  offline_time = []
  N_log = []
  Block_size_log = []

  #For each mode we have to get the time for each n in N, for all the block_sizes
  for b in BLOCK_SIZE:   
    bs_exists = False
    log_file_name = RESULTS_FOLDER+"/"+str(m)+"_"+str(b)+".lg"
    log_file_exists = os.path.isfile(log_file_name)
    if(log_file_exists):
      log_file = open(log_file_name, "r")
      bs_exists = True
    else:
      continue

    for line in log_file:
      values = line.split(',')
      n = int(values[0].strip()) 
      if(n in N):
        if(bs_exists):
          N_log.append(n)
          Block_size_log.append(b)
        ttime.append(float(values[1].strip()))
        ttime_stddev.append(float(values[2].strip()))
        pt = float(values[3].strip())/time_divisor 
        ptime.append(pt)
        ptime_stddev.append(float(values[4].strip())/time_divisor)
        # Last layer time (applicapble to BORPStream (3 and 13)
        # For Waksman Network modes this corresponds to generate permutation time
        ll = float(values[5].strip())/time_divisor 
        lltime.append(ll)
        lltime_stddev.append(float(values[6].strip())/time_divisor)
        if(m in {7,9,10,15,17}):
          gp = float(values[5].strip())/time_divisor
          time_gp.append(gp)
          cb = float(values[8].strip())/time_divisor 
          time_cb.append(cb)
          ap = float(values[10].strip())/time_divisor
          time_ap.append(ap)
          OSWAP_count_gp.append(int(values[12].strip()))
          OSWAP_count_cb.append(int(values[13].strip()))
          OSWAP_count_ap.append(int(values[14].strip()))
          if(m in {7,9}):
            offline_time.append(gp + cb)
            online_time.append(ap)
        #OSWAP_count = OSWAP_count_ap for MODES 6,7,9
        if(m in {7,9,10,15,17}):
          OSWAP_count.append(int(values[12].strip()) + int(values[13].strip()) + int(values[14].strip()))
        else:
          OSWAP_count.append(int(values[7].strip()))

        if(m == 3 or m ==13):
          # For BORP modes:
	  # OSWAP_count_cb = Phase1 OSWAP count 
          # OSWAP_count_ap = Phase2 OSWAP count 
          OSWAP_count_cb.append(int(values[8].strip()))
          OSWAP_count_ap.append(int(values[9].strip()))
        if(m==3):
          online_time.append(ll)
          offline_time.append(pt - ll)
        if(m == 13):
          qt = float(values[10].strip())/time_divisor
          qsort_time.append(qt)
          qsort_time_stddev.append(float(values[11].strip())/time_divisor)
          online_time.append(ll + qt)
          offline_time.append(pt - (ll + qt))
        if(m == 14):
          qsort_time.append(float(values[8].strip())/time_divisor)
          qsort_time_stddev.append(float(values[9].strip())/time_divisor)
        if(m == 15):
          qt = float(values[15].strip())/time_divisor 
          qsort_time.append(qt)
          qsort_time_stddev.append(float(values[16].strip())/time_divisor)
          online_time.append(ap + qt)
          offline_time.append(gp + cb)

  if(PLOT_MODE in {1,3,4,6}):
    X_AXIS = N_log
  else:
    X_AXIS = Block_size_log

  #print(m, label_ctr)
  label_str1 = LABEL_STRINGS[label_ctr]

  if(PLOT_MODE in TIME_PLOT_MODES):
    if(m in ONLY_PTIME_MODES):
      points1 = plt.plot(X_AXIS, ptime, POINT_STYLES[CMT[m]], label=label_str1, color = str("C")+str(CMT[m]), linestyle='solid', linewidth=2.5, mew=4)
    else:
      if(m in {17, 10, 12, 14}):
        points1 = plt.plot(X_AXIS, ptime, POINT_STYLES[CMT[m]], label=label_str1, color = str("C")+str(CMT[m]), linestyle='solid', linewidth=2.5, mew=4)
      else:
        points1 = plt.plot(X_AXIS, ptime, POINT_STYLES[CMT[m]], label=label_str1+": Total", color = str("C")+str(CMT[m]), linestyle='solid', linewidth=2.5, mew=4)
    errbar1 = plt.errorbar(X_AXIS, ptime, ptime_stddev, None,'r', ls='none', mew=3)
  else:
    if(m in ONLY_PTIME_MODES):
      points1 = plt.plot(X_AXIS, OSWAP_count, POINT_STYLES[CMT[m]], label=label_str1, color = str("C")+str(CMT[m]), linestyle='solid', linewidth=2.5, mew=4)
    else:
      points1 = plt.plot(X_AXIS, OSWAP_count, POINT_STYLES[CMT[m]], label=label_str1+": Total", color = str("C")+str(CMT[m]), linestyle='solid', linewidth=2.5, mew=4)


  if(m in QSORT_MODES):
    if(PLOT_QSORT_TIME):
      points1 = plt.plot(X_AXIS, qsort_time, POINT_STYLES[ps_ctr % PS_SIZE], label=label_str1+": QSort", color = str("C")+str(color_ctr%10), linestyle=LINE_STYLES[ls_ctr%len(LINE_STYLES)], linewidth=2.5, mew=4)
      errbar1 = plt.errorbar(X_AXIS, qsort_time, qsort_time_stddev, None,'r', ls='none', mew=3)

  if(m in WAKSMAN_MODES):
    if(PLOT_GEN_PERM):
      points2 = plt.plot(X_AXIS, time_gp, POINT_STYLES[ps_ctr % PS_SIZE], label=label_str1 + ": Gen Perm", color = str("C")+str(color_ctr%10), linestyle=LINE_STYLES[ls_ctr%len(LINE_STYLES)], linewidth=2.5, mew=4)
      errbar2 = plt.errorbar(X_AXIS, ptime, ptime_stddev[str(n)], None,'r', ls='none', mew=3)
    
    if(PLOT_CONTROL_BITS):
      points3 = plt.plot(X_AXIS, time_cb, POINT_STYLES[ps_ctr % PS_SIZE], label=label_str1 + ": Control Bits", color = str("C")+str(color_ctr%10), linestyle=LINE_STYLES[ls_ctr%len(LINE_STYLES)], linewidth=2.5, mew=4)
      errbar3 = plt.errorbar(X_AXIS, ptime, ptime_stddev, None,'r', ls='none', mew=3)

    if(PLOT_APPLY_PERM):
      points4 = plt.plot(X_AXIS, time_ap, POINT_STYLES[ps_ctr % PS_SIZE], label=label_str1 + ": Apply Perm", color = str("C")+str(color_ctr%10), linestyle=LINE_STYLES[ls_ctr%len(LINE_STYLES)], linewidth=2.5, mew=4)
      errbar4 = plt.errorbar(X_AXIS, ptime, ptime_stddev, None,'r', ls='none', mew=3)

  if(m in WAKSMAN_SHUFFLE_MODES or m in BORP_MODES):
    if(PLOT_ONLINE_TIME):
      local_label = ""
      if(m == 7):
        local_label = "WaksOn"
      elif(m == 15):
        local_label = "WaksOn + QS"
      else:
        local_label = label_str1
      local_label=label_str1 + ": Online"
      points5 = plt.plot(X_AXIS, online_time, POINT_STYLES[CMO[m]], label=local_label, color = str("C")+str(CMO[m]), linestyle='dashed', linewidth=2.5, mew=4)

    if(PLOT_OFFLINE_TIME):
      local_label = ""
      if(m in 7 or 15):
        local_label = "WaksOff"
      else:
        label=label_str1 + ": Offline"
      points6 = plt.plot(X_AXIS, offline_time, POINT_STYLES[ps_ctr % PS_SIZE], label=local_label, color = str("C")+str(color_ctr%10), linestyle=LINE_STYLES[ls_ctr%len(LINE_STYLES)], linewidth=2.5, mew=4)


  label_ctr+=1
  ls_ctr+=1

ax = plt.subplot()
ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
#handles, labels = ax.get_legend_handles_labels()

if (X_AXIS_LOGSCALE):
  plt.yscale('log')
if (Y_AXIS_LOGSCALE):
  plt.xscale('log')

# To Convert x-axis to powers of 2 labels


if (PLOT_MODE in {1,4}):
  if(PLOT_MODE==1):
    plt.title('Comparing Shuffling Algorithms\nComputation time (in s) vs number of items', fontdict = {'fontsize' : title_font_size} )
  elif(PLOT_MODE==4):
    plt.title('Comparing Sorting Algorithms\nComputation time (in s) vs number of items', fontdict = {'fontsize' : title_font_size} )

  plt.xlabel('Number of items', fontsize = axis_label_font_size)
  plt.ylabel('Computation time (in s)', fontsize = axis_label_font_size)
  ax.set_xticks(N)
  labels = convertTo2PowLabels(N)
  ax.set_xticklabels(labels)
elif (PLOT_MODE in {2, 5}):
  if(PLOT_MODE==2):
    plt.title('Comparing Shuffling Algorithms\nComputation time (in s) vs item size in bytes', fontdict = {'fontsize' : title_font_size} )
  elif(PLOT_MODE==5):
    plt.title('Comparing Sorting Algorithms\nComputation time (in s) vs item size in bytes', fontdict = {'fontsize' : title_font_size} )
  plt.xlabel('Item size in bytes', fontsize = axis_label_font_size)
  plt.ylabel('Computation time (in s)', fontsize = axis_label_font_size)
  ax.set_xticks(LABELS_Y)
  labels = convertTo2PowLabels(LABELS_Y)
  ax.set_xticklabels(labels)
elif (PLOT_MODE in {3, 6}):
  plt.title('Number of OSWAPs vs number of items', fontdict = {'fontsize' : title_font_size} )
  plt.xlabel('Number of items', fontsize = axis_label_font_size)
  plt.ylabel('Number of OSWAPs', fontsize = axis_label_font_size)
  ax.set_xticks(N)
  labels = convertTo2PowLabels(N)
  ax.set_xticklabels(labels)

ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.25), ncol=2, shadow=True, prop={'size': legend_font_size})

plt.minorticks_off()
plt.savefig(save_location, bbox_inches='tight')
plt.close()
