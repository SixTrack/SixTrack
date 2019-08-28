#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting beam overlap between beam 1 and beam 2
Veronica Berglyd Olsen, BE-ABP-HSS
"""

import numpy as np
import matplotlib.pyplot as plt
from os import path

elName = ""
iTurn  = 0
xLim   = [0,0]
yLim   = [0,0]
xNum   = 0
yNum   = 0

# Collect data from PDF file
iLine = 0
nPart = 0
with open("scatter_pdf.dat",mode="r") as inFile:
  for inLine in inFile:
    iLine += 1
    if iLine == 1:
      spLine  = inLine.split()
      elName  = spLine[-1]
    elif iLine == 2:
      spLine  = inLine.split()
      iTurn   = int(spLine[-1])
    elif iLine == 3:
      spLine  = inLine.split()
      xLim[0] = float(spLine[-3])
      xLim[1] = float(spLine[-2])
      xNum    = int(spLine[-1])
    elif iLine == 4:
      spLine  = inLine.split()
      yLim[0] = float(spLine[-3])
      yLim[1] = float(spLine[-2])
      yNum    = int(spLine[-1])

bPDF = np.zeros((xNum, yNum))
xPos = 0
yPos = 0
with open("scatter_pdf.dat",mode="r") as inFile:
  for inLine in inFile:
    iLine += 1
    if inLine[:1] == "#":
      continue
    theVals = inLine.split()
    for aVal in theVals:
      bPDF[yPos,xPos] = float(aVal)
      xPos += 1
      if xPos >= xNum:
        xPos = 0
        yPos += 1
      if yPos > yNum:
        print("ERROR PDF map index out of bounds")
        exit(1)

print("Element: %s" % elName)
print("Turn:    %d" % iTurn)
print("X-Range: %13.6f %13.6f %d" % (xLim[0], xLim[1], xNum))
print("Y-Range: %13.6f %13.6f %d" % (yLim[0], yLim[1], yNum))

iLine = 0
nPart = 0
iPart = 0
bScat = None
with open("scatter_density.dat",mode="r") as inFile:
  for inLine in inFile:
    iLine  += 1
    theVals = inLine.split()
    if iLine == 3:
      nPart = int(theVals[-1])
      bScat = np.zeros((nPart,2))
      print("NPart:   %d" % nPart)
      continue
    if inLine[:1] == "#":
      continue
    if len(theVals) != 7:
      print("WARNING Line %d had %d values, expected 7" % (iLine, len(theVals)))
      continue
    bScat[iPart,0] = float(theVals[4])
    bScat[iPart,1] = float(theVals[5])
    iPart += 1

theFig = plt.figure(figsize=(7, 6),dpi=100)
theFig.clf()

plt.imshow(bPDF,extent=[xLim[0],xLim[1],yLim[0],yLim[1]],aspect="auto",zorder=1)
plt.scatter(bScat[:,0],bScat[:,1],1,"#000000",marker=",",zorder=2)
plt.title("Beam 2 at '%s' on turn %d" % (elName, iTurn))
plt.xlim((xLim[0],xLim[1]))
plt.ylim((yLim[0],yLim[1]))
plt.show(block=True)
