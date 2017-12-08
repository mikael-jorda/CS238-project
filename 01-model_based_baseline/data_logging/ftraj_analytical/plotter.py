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

# plt.figure(1)
plt.figure(num=1, figsize=(8, 8), dpi=80)
plt.subplot(2,1,1)
plt.plot(rtime, rforce[:,1], label="Actual Force (z)", color="k", linewidth="1.5")
plt.plot(etime, eforce[:,1], label="Estimated Force (z)", color="r", linewidth="1.5", linestyle="--")
plt.plot(rtime, rforce[:,0], label="Actual Force (y)", color="g", linewidth="1.5")
plt.plot(etime, eforce[:,0], label="Estimated Force (z)", color="b", linewidth="1.5", linestyle="--")
plt.legend(bbox_to_anchor=(1.0, 0.8), loc=1, borderaxespad=0.5)
plt.title("Hold Position, Analytical approach\nForce", fontsize=20)
plt.ylabel("Force (N)", fontsize=17)
# plt.legend()


# plt.figure(2)
plt.subplot(2,1,2)
plt.plot(rtime, rpos[:,1], label="Actual contact position (z)", color="k", linewidth="1.5")
plt.plot(etime, epos[:,1], label="Estimated contact position (z)", color="r", linewidth="1.5", linestyle="--")
plt.plot(rtime, rpos[:,0], label="Actual contact position (y)", color="g", linewidth="1.5")
plt.plot(etime, epos[:,0], label="Estimated contact position (y)", color="b", linewidth="1.5", linestyle="--")
plt.legend(bbox_to_anchor=(1.0, 0.5), loc=1, borderaxespad=0.5)
plt.title("Contact Position", fontsize=20)
plt.xlabel("time (s)", fontsize=17)
plt.ylabel("Position (m)", fontsize=17)
# plt.legend()



plt.figure(3)
plt.plot(ptime, pdesired[:,0], label="yd")
plt.plot(ptime, pdesired[:,1], label="zd")
plt.plot(ptime, pactual[:,0], label="y")
plt.plot(ptime, pactual[:,1], label="z")
plt.legend()


plt.show()


