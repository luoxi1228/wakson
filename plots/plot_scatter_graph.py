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
# PLOT_MODE = 1 for scatter plot of shuffle algorithms
# PLOT_MODE = 2 for scatter plot of sort algorithms
PLOT_MODE = 2

RESULTS_FOLDER = "Our_Data"

# Set scale of axis:
X_AXIS_LOGSCALE = True
Y_AXIS_LOGSCALE = True

PLOT_TIME_IN_MS = False

BLOCK_SIZE = [256]
#BLOCK_SIZE = [4096]

###############################################################################

if(len(sys.argv)>1):
  RESULTS_FOLDER = (sys.argv[1])
  PLOT_MODE = int(sys.argv[2])
  bs_mode = int(sys.argv[3])
  if(bs_mode==1):
    BLOCK_SIZE = [256]
  elif(bs_mode==2):
    BLOCK_SIZE = [4096]

###############################################################################

# Select other parameters depending on CONFIGS
flag_BV1 = True
if(not(os.path.isfile(RESULTS_FOLDER+"/"+"33_256.lg"))):
  flag_BV1 = False

if (PLOT_MODE == 1):
  MODE = [5, 2, 3, 7, 9, 33]
  N = [1048576]
  LABEL_STRINGS = ["BitonicShuffle", "ORShuffle", "BORPStream (V2)", "WaksShuffle", "Nassimi-Sahni", "BORPStream (V1)"]
  if(not(flag_BV1)):
    MODE = [5, 2, 3, 7, 9]
    LABEL_STRINGS = ["BitonicShuffle", "ORShuffle", "BORPStream (V2)", "WaksShuffle", "Nassimi-Sahni"]
  bs_val = str(BLOCK_SIZE[0])
  save_location = RESULTS_FOLDER+ '/Scatter_Shuffles_'+bs_val+'.png'

elif (PLOT_MODE == 2):
  MODE = [10, 12, 13, 14, 15, 17, 43]
  if(BLOCK_SIZE[0]==256):
    LABEL_STRINGS = ["Nassimi-Sahni", "ORShuffle\n + QS", "BORPStream (V2) + QS", "Bitonic Sort", "WaksShuffle + QS", "WaksSort",  "BORPStream (V1) + QS"] 
  else:
    LABEL_STRINGS = ["Nassimi-Sahni", "ORShuffle\n + QS", "BORPStream\n (V2) + QS", "Bitonic Sort", "WaksShuffle + QS", "WaksSort",  "BORPStream\n (V1) + QS"]
  if(not(flag_BV1)):
    MODE = [10, 12, 13, 14, 15, 17]
    if(BLOCK_SIZE[0]==256):
      LABEL_STRINGS = ["Nassimi-Sahni", "ORShuffle\n + QS", "BORPStream (V2) + QS", "Bitonic Sort", "WaksShuffle + QS", "WaksSort"] 
    else:
      LABEL_STRINGS = ["Nassimi-Sahni", "ORShuffle\n + QS", "BORPStream\n (V2) + QS", "Bitonic Sort", "WaksShuffle + QS", "WaksSort"]

  N = [1048576]
  bs_val = str(BLOCK_SIZE[0])
  save_location = RESULTS_FOLDER+ '/Scatter_Sorts_'+bs_val+'.png'

time_divisor = 1000.0
if(PLOT_TIME_IN_MS):
  time_divisor = 1.0
  
###############################################################################

POINT_STYLES = [".","v","s","+","^","d",">","x","*","o","D"]
LINE_STYLES = ['solid', 'dotted', 'dashed', 'dashdot']

PS_SIZE = len(POINT_STYLES)
LS_SIZE = len(LINE_STYLES)
# Colors range from C0 to C9

#OG:
#COLOR_MAP_TOTAL 
CMT = {5:8, 2:1, 3:4, 7:5, 9:6, 10:6, 12:8, 13:4, 14:1, 15:5, 17:9, 33:0, 43:0}
#COLOR_MAP_ONLINE 
CMO = {3:4, 7:2, 9:6, 15:2, 13:4}

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
title_font_size = 20
legend_font_size = 18
tick_font_size = 18
axis_label_font_size = 20

###############################################################################

TIME_PLOT_MODES = {1, 2, 4, 5}
OSWAP_PLOT_MODES = {3, 6}

ONLY_PTIME_MODES = [1, 2, 4, 5, 10, 12, 14, 17, 33, 43]
BORP_MODES = [3, 13]
WAKSMAN_MODES = [7, 9, 10, 15, 17]
WAKSMAN_SHUFFLE_MODES = [7, 9, 15]
WAKSMAN_SORT_MODE = [10, 17]
QSORT_MODES = [13, 14, 15]

###############################################################################

label_ctr=0
color_ctr=0
ps_ctr=0
ls_ctr=0


plt.figure(figsize=(8.5,3))
plt.grid(True)

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

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

  label= LABEL_STRINGS[label_ctr]

  if(PLOT_MODE == 2 and BLOCK_SIZE[0] == 256):
    if(m not in ONLY_PTIME_MODES):
      plt.scatter(ptime[0], online_time[0], marker = POINT_STYLES[CMO[m]], color = str("C")+str(CMO[m]), s = 175.0)
      if(m==13):
        plt.text(ptime[0] + 1, online_time[0] - 0.2, label, fontsize = legend_font_size)
      elif(m==15):
        plt.text(ptime[0] + 2, online_time[0] - 0.1, label, fontsize = legend_font_size)
    else:
      plt.scatter(ptime[0], ptime[0], marker = POINT_STYLES[CMT[m]], color = str("C")+str(CMT[m]), s = 175.0)
      if(m==10):
        plt.text(ptime[0] - 200, ptime[0] - 50, label, fontsize = legend_font_size)
      elif(m==12):
        plt.text(ptime[0], ptime[0] + 1.5, label, fontsize = legend_font_size, rotation = 'vertical')
      elif(m==14):
        plt.text(ptime[0] - 0.25, ptime[0] - 1.2, label, fontsize = legend_font_size)
      elif(m==17):
        plt.text(ptime[0] + 3, ptime[0] -1.5, label, fontsize = legend_font_size)
      elif(m==43):
        plt.text(ptime[0] + 1, ptime[0] - 2.5, label, fontsize = legend_font_size)

  elif(PLOT_MODE == 2 and BLOCK_SIZE[0] == 4096):
    if(m not in ONLY_PTIME_MODES):
      plt.scatter(ptime[0], online_time[0], marker = POINT_STYLES[CMO[m]], color = str("C")+str(CMO[m]), s = 175.0)
      if(m==13):
        plt.text(ptime[0] + 4, online_time[0] - 3, label, fontsize = legend_font_size)
      elif(m==15):
        plt.text(ptime[0] + 2, online_time[0] - 0.6, label, fontsize = legend_font_size)
    else:
      plt.scatter(ptime[0], ptime[0], marker = POINT_STYLES[CMT[m]], color = str("C")+str(CMT[m]), s = 175.0)
      if(m==10):
        plt.text(ptime[0] - 155, ptime[0] - 40, label, fontsize = legend_font_size)
      elif(m==12):
        plt.text(ptime[0], ptime[0] + 15, label, fontsize = legend_font_size, rotation = 'vertical')
      elif(m==14):
        plt.text(ptime[0] -2, ptime[0] - 11, label, fontsize = legend_font_size)
      elif(m==17):
        plt.text(ptime[0] + 2, ptime[0] - 1.5, label, fontsize = legend_font_size)
      elif(m==43):
        plt.text(ptime[0] + 10, ptime[0] - 25, label, fontsize = legend_font_size)

  elif(PLOT_MODE == 1 and BLOCK_SIZE[0] == 256):
    if(m not in ONLY_PTIME_MODES):
      plt.scatter(ptime[0], online_time[0], marker = POINT_STYLES[CMO[m]], color = str("C")+str(CMO[m]), s = 175.0)
      if(m==3):
        plt.text(ptime[0] + 1, online_time[0] - 0.05, label, fontsize = legend_font_size)
      elif(m==7):
        plt.text(ptime[0] + 1, online_time[0] -0.03, label, fontsize = legend_font_size)
      elif(m==9):
        plt.text(ptime[0] - 200, online_time[0] - 0.1, label, fontsize = legend_font_size)
    else:
      plt.scatter(ptime[0], ptime[0], marker = POINT_STYLES[CMT[m]], color = str("C")+str(CMT[m]), s = 175.0)
      if(m==2):
        plt.text(ptime[0] + 0.5, ptime[0]-0.4, label, fontsize = legend_font_size)
      elif(m==5):
        plt.text(ptime[0] + 0.5, ptime[0], label, fontsize = legend_font_size)
      elif(m==33):
        plt.text(ptime[0] + 1, ptime[0] -0.5, label, fontsize = legend_font_size)

  elif(PLOT_MODE == 1 and BLOCK_SIZE[0] == 4096):
    if(m not in ONLY_PTIME_MODES):
      plt.scatter(ptime[0], online_time[0], marker = POINT_STYLES[CMO[m]], color = str("C")+str(CMO[m]), s = 175.0)
      if(m==3):
        plt.text(ptime[0] - 87, online_time[0]- 1 , label, fontsize = legend_font_size)
      elif(m==7):
        plt.text(ptime[0] + 1, online_time[0] - 0.6, label, fontsize = legend_font_size)
      elif(m==9):
        plt.text(ptime[0] - 153, online_time[0] - 0.5, label, fontsize = legend_font_size)
    else:
      plt.scatter(ptime[0], ptime[0], marker = POINT_STYLES[CMT[m]], color = str("C")+str(CMT[m]), s = 175.0)
      if(m==2):
        plt.text(ptime[0] + 0.5, ptime[0] - 9.5, label, fontsize = legend_font_size)
      elif(m==5):
        plt.text(ptime[0] + 3, ptime[0] - 1.5, label, fontsize = legend_font_size)
      elif(m==33):
        plt.text(ptime[0] - 72, ptime[0] - 5, label, fontsize = legend_font_size)

  label_ctr+=1

if (X_AXIS_LOGSCALE):
  plt.xscale('log')
if (Y_AXIS_LOGSCALE):
  plt.yscale('log')

ax = plt.subplot()
ax.tick_params(axis='both', which='major', labelsize=tick_font_size)

if(PLOT_MODE == 1 and BLOCK_SIZE[0] == 256):
  x_ticks = [2, 4, 8, 16, 32, 64, 128, 256]
  x_labels = [r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  y_ticks = [0.5, 1, 2, 4, 8,]
  y_labels = [r'$2^{-1}$', r'$2^0$',r'$2^1$', r'$2^2$', r'$2^3$']
  ax.set_xticks(x_ticks)
  ax.set_yticks(y_ticks)
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)
  plt.ylim(0.4,9)

elif(PLOT_MODE==1 and BLOCK_SIZE[0] == 4096):
  x_ticks = [2, 4, 8, 16, 32, 64, 128, 256]
  x_labels = [r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  y_ticks = [1, 2, 4, 8, 16, 32, 64, 128, 256]
  y_labels = [ r'$2^0$', r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  ax.set_xticks(x_ticks)
  ax.set_yticks(y_ticks)
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)
  plt.ylim(6,140)
  plt.xlim(15, 280)

elif(PLOT_MODE == 2 and BLOCK_SIZE[0] == 256):
  x_ticks = [2, 4, 8, 16, 32, 64, 128, 256]
  x_labels = [r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  y_ticks = [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
  y_labels = [r'$2^{-1}$', r'$2^0$',r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  ax.set_xticks(x_ticks)
  ax.set_yticks(y_ticks)
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)

elif(PLOT_MODE == 2 and BLOCK_SIZE[0] == 4096):
  #x_labels = [2.5, 5, 10, 20, 40, 80, 160]
  x_ticks = [2, 4, 8, 16, 32, 64, 128, 256]
  #x_labels = [2, 4, 8, 16, 32, 64, 128, 256]
  x_labels = [r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  y_ticks = [1, 2, 4, 8, 16, 32, 64, 128, 256]
  y_labels = [ r'$2^0$', r'$2^1$', r'$2^2$', r'$2^3$', r'$2^4$', r'$2^5$', r'$2^6$', r'$2^7$', r'$2^8$']
  ax.set_xticks(x_ticks)
  ax.set_yticks(y_ticks)
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)
  #plt.ylim(6,140)
  plt.xlim(15, 280)



if(PLOT_TIME_IN_MS):
  plt.title('Sorting: Online and Total Time', fontdict = {'fontsize' : title_font_size} )
  plt.xlabel('Total Time (in ms)', fontsize = axis_label_font_size)
  plt.ylabel('Online Time (in ms)', fontsize = axis_label_font_size)
else:
  if(PLOT_MODE==1):
    if(BLOCK_SIZE[0] == 4096):
      plt.title('Shuffling: Online and Total Time (%s = 4 KiB)' % ('\u03B2'), fontdict = {'fontsize' : title_font_size} )
    elif(BLOCK_SIZE[0] == 256):
      plt.title('Shuffling: Online and Total Time (%s = 256 B)' % ('\u03B2'), fontdict = {'fontsize' : title_font_size} )
  elif(PLOT_MODE==2):
    if(BLOCK_SIZE[0] == 4096):
      plt.title('Sorting: Online and Total Time (%s = 4 KiB)' % ('\u03B2'), fontdict = {'fontsize' : title_font_size} )
    elif(BLOCK_SIZE[0] == 256):
      plt.title('Sorting: Online and Total Time (%s = 256 B)'% ('\u03B2'), fontdict = {'fontsize' : title_font_size} )
  plt.xlabel('Total Time (in s)', fontsize = axis_label_font_size)
  plt.ylabel('Online Time (in s)', fontsize = axis_label_font_size)


#plt.xlim(1,500000)
#plt.ylim(1,500000)
add_identity(ax, color='0.8', ls='--')
plt.minorticks_off()
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=3, shadow=True, prop={'size': legend_font_size})

plt.savefig(save_location, bbox_inches='tight')
plt.close()

