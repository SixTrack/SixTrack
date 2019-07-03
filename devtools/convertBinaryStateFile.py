#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SixTrack DevTools

  SixTrack DevTools :: Convert Binary Particle State File
 =========================================================
  Converts a binary particle state file to a text file
  By: Veronica Berglyd Olsen
      CERN (BE-ABP-HSS)
      Geneva, Switzerland

"""

import sys
from os import path
from struct import unpack

if len(sys.argv) != 2:
  print("\n"
    "convertBinaryStateFile.py takes one input argument:\n"
    " - filename : The file to plot. Must be initial_state.dat or final_state.dat (no binary files).\n"
  )
  exit(0)

inFile = sys.argv[1]

if not path.isfile(inFile):
  raise FileNotFoundError("No file named '%s'" % inFile)

inData  = open(inFile, "rb")
outData = open(inFile+".dat", "w")

napxo,napx,npart,numl = unpack("<iiii", inData.read(16))

outData.write("# Tracking\n")
outData.write("# NPart Start     = %d\n" % napxo)
outData.write("# NPart End       = %d\n" % napx)
outData.write("# NPart Allocated = %d\n" % npart)
outData.write("# NTurns          = %d\n" % numl)

nucm0,e0,e0f,aa0,zz0,qq0,iions = unpack("<dddhhhh", inData.read(32))

outData.write("#\n")
outData.write("# Reference Particle\n")
outData.write("# Mass [MeV]      = %23.16e\n" % nucm0)
outData.write("# Energy [MeV]    = %23.16e\n" % e0)
outData.write("# Momentum [MeV]  = %23.16e\n" % e0f)
outData.write("# Atomic Mass     = %d\n"      % aa0)
outData.write("# Atomic Number   = %d\n"      % zz0)
outData.write("# Charge          = %d\n"      % qq0)

tmp = [0.0]*8

outData.write("#\n")
outData.write("# Closed Orbit [x, xp, y, yp, sigma, dp]\n")
tmp[:4] = unpack("<dddd", inData.read(32))
outData.write("# 4D Closed Orbit = %s\n" % (" ".join(f"{n:23.16e}" for n in tmp[:4])))
tmp[:6] = unpack("<dddddd", inData.read(48))
outData.write("# 6D Closed Orbit = %s\n" % (" ".join(f"{n:23.16e}" for n in tmp[:6])))

outData.write("#\n")
tmp[:3] = unpack("<ddd", inData.read(24))
outData.write("# Tune            = %s\n" % (" ".join(f"{n:23.16e}" for n in tmp[:3])))

for i in range(1,7):
  tmp[:6] = unpack("<dddddd", inData.read(48))
  outData.write("# TAS(%d,1:6)      = %s\n" % (i," ".join(f"{n:23.16e}" for n in tmp[:6])))

outData.write("#\n")
outData.write("# partID parentID lost prim                       x                       y                      xp                      yp                   sigma                      dp                       p                       e")
if iions == 1:
  outData.write("                    mass    A    Z    Q       PDGid\n")
else:
  outData.write("\n")


for i in range(npart):
  partID,parentID,iLost,iPrim = unpack("<iiii", inData.read(16))
  outData.write("%8d %8d" % (partID,parentID))
  if iLost == 1:
    outData.write("    T")
  else:
    outData.write("    F")
  if iPrim == 1:
    outData.write("    T")
  else:
    outData.write("    F")
  tmp[:8] = unpack("<dddddddd", inData.read(64))
  outData.write(" %s" % (" ".join(f"{n:23.16e}" for n in tmp[:8])))
  if iions == 1:
    nucm,naa,nzz,nqq,pdgid,dum1,dum2 = unpack("<dhhhihi", inData.read(24))
    outData.write(" %23.16e %4d %4d %4d %11d" % (nucm,naa,nzz,nqq,pdgid))

  outData.write("\n")

inData.close()
outData.close()
