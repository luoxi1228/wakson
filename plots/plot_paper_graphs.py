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

RESULTS_FOLDER = "Our_Data"
if(len(sys.argv)>2):
  print("Incorrect usage of plotter script.\n Expected ./plot_paper_graphs.py <Experiment_results_folder>")
elif(len(sys.argv)==2):
  RESULTS_FOLDER = sys.argv[1]	

###############################################################################

subprocess.run(["python3", "plot_graph.py", RESULTS_FOLDER, "2"])
subprocess.run(["python3", "plot_graph.py", RESULTS_FOLDER, "5"])

subprocess.run(["python3", "plot_scatter_graph.py", RESULTS_FOLDER, "1", "1"])
subprocess.run(["python3", "plot_scatter_graph.py", RESULTS_FOLDER, "1", "2"])
subprocess.run(["python3", "plot_scatter_graph.py", RESULTS_FOLDER, "2", "1"])
subprocess.run(["python3", "plot_scatter_graph.py", RESULTS_FOLDER, "2", "2"])

###############################################################################
