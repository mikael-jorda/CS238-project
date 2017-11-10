#!/usr/bin/env python

# read text files and plot them

import matplotlib.pyplot as plt
import numpy as np
import sys

# data file to read given as argument
if len(sys.argv) < 2:
	print "Give the name of the file to read as an argument\n"
	exit()

file = np.loadtxt(sys.argv[1] ,skiprows=1)

time = file[:,0]
commanded_torques = file[:,1:5]
# cmap = ['r','k']
# 

norm1_torques = abs(commanded_torques)
total_power = np.sum(commanded_torques,1)

plt.figure(1)
plt.plot(commanded_torques,label="fgc")
plt.legend()

plt.figure(2)
plt.plot(total_power)
plt.title("Total motor power")


plt.show()


