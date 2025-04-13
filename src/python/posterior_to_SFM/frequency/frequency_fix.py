import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

freqs = np.loadtxt('./frequency.csv', dtype = np.float64, delimiter=' ', unpack=True)

out_path = '/home/jeremy/Documents/GRSW/frequency_fixed.csv'
out_file = open(out_path, "w")

N = freqs.shape[1]

for i in range(10500):
      row           = freqs[:,i]
      before        = row[1399]
      after         = row[1401]
      freqs[1400,i] = 0.5*(before + after)

for i in range(N):
      row = freqs[:,i]
      n   = len(row)
      for j in range(n - 1):
            out_file.write(str(row[j]) + ",")
      out_file.write(str(row[-1]) + "\n")
      print(i,'/',N)

out_file.close()











