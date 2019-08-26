#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SixTrack DevTools

  SixTrack DevTools :: Plot Particle State
 ==========================================
  Script for plotting the content of initial_state.dat or final_state.dat files
  By: Veronica Berglyd Olsen
      CERN (BE-ABP-HSS)
      Geneva, Switzerland

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from os import path
from math import ceil

if len(sys.argv) < 3:
  print("\n"
    "plotParticleState.py takes two or three input arguments:\n"
    " - filename : The file to plot. Must be initial_state.dat or final_state.dat (no binary files).\n"
    " - column   : Comma separated list of columns to plot (no spaces between them).\n"
    "              Available: x,y,xp,yp,sigma,dp,p,e\n"
    " - saveto   : [Optional] file name to save figure to" 
  )
  exit(0)

# Constants

colMap = {"x":4,"y":5,"xp":6,"yp":7,"sigma":8,"dp":9,"p":10,"e":11}

# Save the arguments and check that they make sense

thePath = sys.argv[1]

if not path.isfile(thePath):
  print("File not found: %s" % thePath)
  exit(1)

theFile = path.basename(thePath)
theCols = []
theData = {}

for aCol in sys.argv[2].split(","):
  aCol = aCol.strip() # Just to make sure!
  if aCol in colMap.keys():
    theCols.append(aCol)
    theData[aCol] = []
  else:
    print("Not a valid column: %s" % aCol)
    exit(1)

if len(sys.argv) > 3:
  saveTo = sys.argv[3]
else:
  saveTo = None

# Write a summary
print((
  "\n"
  "Will plot the following:\n"
  " - File:    {filename:s}\n"
  " - Columns: {columns:s}\n"
).format(
  filename = thePath,
  columns  = ", ".join(theCols)
))

# Collect data from file
iLine = 0
nPart = 0
with open(thePath,mode="r") as inFile:
  for inLine in inFile:
    iLine += 1
    if inLine[:1] == "#":
      continue
    theVals = inLine.split()
    if len(theVals) < 12:
      print("WARNING: Line %d has fewer columns than expected. Skipping." % iLine)
      continue
    nPart += 1
    for aCol in theCols:
      try:
        theData[aCol].append(float(theVals[colMap[aCol]]))
      except:
        print("ERROR: Could not convert to float value on line %d, column %d" % (iLine,colMap[aCol]))

print("Read %d lines from file." % iLine)
print("Read %d particles from file." % nPart)

# Plot

iFig   = 0
nBins  = 100
theFig = plt.figure(iFig,figsize=(16, 9),dpi=100)
nRow = 2
nCol = ceil(len(theCols)/2)
theFig.clf()
for aCol in theCols:
  iFig  += 1
  dHist, binEdges   = np.histogram(theData[aCol], bins=nBins, density=True)
  binCentres        = (binEdges[:-1] + binEdges[1:])/2
  
  sFig1 = plt.subplot(nRow,nCol,iFig)
  plt.step(binCentres,dHist,where="mid")
  plt.title("Distribution of %s" % aCol)

plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
if saveTo is not None:
  plt.savefig(saveTo)
plt.show(block=True)


