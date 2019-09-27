#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) == 2:
  fName = sys.argv[1]
else:
  print("ERROR plotRand.py takes a file name as input parameter.")
  exit(1)

if not os.path.isfile(fName):
  raise FileNotFoundError("Missing file: %s" % fName)

aData = []
with open(fName, mode="r") as inFile:
  for inLine in inFile:
    aData.append(float(inLine))

aData = np.array(aData)
nSamp = len(aData)

# Plot
figOne = plt.figure(1,figsize=(10, 8),dpi=100)
figOne.clf()

pHist, pG = np.histogram(aData, bins=200, density=True)
pC = (pG[:-1] + pG[1:])/2
plt.step(pC,pHist,where="mid",label="Mean: %6.3f\nStd: %6.3f" % (np.mean(aData), np.std(aData)))
plt.legend()
plt.title("File: %s" % fName)
plt.savefig(fName[:-4]+".png")
# plt.ylim((0.0,1.5))

plt.show(block=True)
