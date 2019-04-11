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

def runTests(theTests, theArgs, nLines, nSkip):

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

    with open("error_results.log","a") as outFile:
      outFile.write("\n"*4)
      outFile.write(" Test: %s\n" % aTest)
      outFile.write(("#"*132)+"\n")
      outFile.write("\n")

    theBlock = "\n".join(theTests[aTest])
    with open("fort.3","w") as outFile:
      outFile.write(tmpF3.replace("%INIT%",theBlock))
    exCode  = os.system("%s > fort.6" % theArgs[1])
    outBuf += "%-32s %s\n" % (aTest,str(exCode != 0))
    os.system("tail -n%d fort.6 | head -n%d >> error_results.log" % (nLines+nSkip, nLines))

  with open("error_summary.log","w") as outFile:
    outFile.write(outBuf)

  exCode  = os.system("diff -q error_summary.log error_summary.log.canonical")
  exCode += os.system("diff -q error_results.log error_results.log.canonical")

  return int(exCode != 0)
