#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SixTrack Error Test Tools

  SixTrack Error Test Tools
 ===========================
  By: Veronica Berglyd Olsen
      CERN (BE-ABP-HSS)
      Geneva, Switzerland

"""

import os
import subprocess

def runTests(theTests, theArgs):

  if len(theArgs) != 2:
    print("ERROR: Wrong number of arguments.")
    return 1

  try:
    with open("fort.3.template","r") as inFile:
      tmpF3 = inFile.read()
  except:
    print("ERROR: Could not open file fort.3.template")
    return 1

  with open("error_results.log","w") as outFile:
    outFile.write("\n")
    outFile.write(("#"*132)+"\n")
    outFile.write(" SixTrack Error Tests\n")
    outFile.write(("#"*132)+"\n")

  outBuf  = "%-32s %s\n" % ("Test Name","Pass")
  outBuf += ("="*40)+"\n"

  for aTest in sorted(theTests.keys()):
    theBlock = "\n".join(theTests[aTest])
    with open("fort.3","w") as outFile:
      outFile.write(tmpF3.replace("%ERRORTESTS%",theBlock))
    stdOut, stdErr, exCode = sysCall(theArgs[1])
    stdErr = stdErr.replace("\n1\n","\n")       # Removes the final line in stderr for ifort
    stdErr = stdErr.replace("\nSTOP 1\n","\n")  # Removes the final line in stderr for gfortran
    stdErr = stdErr.replace("\nSTOP: 1\n","\n") # Removes the final line in stderr for nagfor
    outBuf += "%-32s %s\n" % (aTest,str(exCode != 0))
    with open("error_results.log","a") as outFile:
      outFile.write("\n"*4)
      outFile.write(" Test: %s\n" % aTest)
      outFile.write(("#"*132)+"\n")
      outFile.write("\n")
      outFile.write(stdErr)

  with open("error_summary.log","w") as outFile:
    outFile.write(outBuf)

  exCode  = os.system("diff -q error_summary.log error_summary.log.canonical")
  exCode += os.system("diff -q error_results.log error_results.log.canonical")

  return int(exCode != 0)

# Wrapper function for system calls
def sysCall(callStr):
  sysP = subprocess.Popen([callStr], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  stdOut, stdErr = sysP.communicate()
  return stdOut.decode("utf-8"), stdErr.decode("utf-8"), sysP.returncode
