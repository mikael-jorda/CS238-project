#!/usr/bin/env python

# read text files and plot them

import matplotlib.pyplot as plt
import numpy as np
import sys

# data files to read
real_file = "real.txt";
estimate_file = "estimates.txt"
position_file = "pos.txt"


rdata = np.loadtxt(real_file ,skiprows=1)
edata = np.loadtxt(estimate_file ,skiprows=1)
pdata = np.loadtxt(position_file, skiprows=1)

rtime = rdata[:,0]
rforce = rdata[:,1:3]
rpos = rdata[:,3::]

etime = edata[:,0]
eforce = edata[:,1:3]
epos = edata[:,3::]

ptime = pdata[:,0]
pdesired = pdata[:,1:3]
pactual = pdata[:,3::]

plt.figure(1)
plt.plot(rtime, rforce[:,0], label="Fry")
plt.plot(rtime, rforce[:,1], label="Frz")
plt.plot(etime, eforce[:,0], label="Fey")
plt.plot(etime, eforce[:,1], label="Fez")
plt.legend()


plt.figure(2)
plt.plot(rtime, rpos[:,0], label="Pry")
plt.plot(rtime, rpos[:,1], label="Prz")
plt.plot(etime, epos[:,0], label="Pey")
plt.plot(etime, epos[:,1], label="Pez")
plt.legend()


plt.figure(3)
plt.plot(ptime, pdesired[:,0], label="yd")
plt.plot(ptime, pdesired[:,1], label="zd")
plt.plot(ptime, pactual[:,0], label="y")
plt.plot(ptime, pactual[:,1], label="z")
plt.legend()


plt.show()


