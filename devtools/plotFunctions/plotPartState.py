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

if len(sys.argv) != 3:
  print("\n"
    "plotParticleState.py takes two input arguments:\n"
    " - filename : The file to plot. Must be initial_state.dat or final_state.dat (no binary files).\n"
    " - column   : Comma separated list of columns to plot (no spaces between them).\n"
    "              Available: x,y,xp,yp,sigma,dp,p,e\n"
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
if not theFile[:15] == "final_state.dat" and not theFile[:17] == "initial_state.dat":
  print("Not a particle state file: %s" % theFile)
  exit(1)

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

iFig  = 0
nBins = 40
for aCol in theCols:
  iFig  += 1
  theFig = plt.figure(iFig,figsize=(6, 4),dpi=100)
  theFig.clf()
  dHist, binEdges   = np.histogram(theData[aCol], bins=nBins, density=True)
  binCentres        = (binEdges[:-1] + binEdges[1:])/2
  # fCoeff, varMatrix = curve_fit(fitGauss, binCentres, dHist, p0=p0)
  # hFit              = fitGauss(binCentres, *fCoeff)
  
  # sFig1 = plt.subplot(2,3,i+1)
  plt.step(binCentres,dHist,where="mid")
  # plt.hist(theData[aCol], bins=nBins, density=True)
  # plt.plot(binCentres,hFit)
  plt.title("Distribution of %s" % aCol)
  # plt.xlabel("mu = %.2f, sigma = %.2f" % (fCoeff[1],fCoeff[2]))

plt.show(block=True)


