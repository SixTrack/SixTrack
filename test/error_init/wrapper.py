#!/usr/bin/env python

import sys, os

defLn1 = (0,0.0,0.0,0.0,0)
defLnN = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,7.0e9,7.0e9,7.0e9)

theTests = {
  "1. Valid INIT"     : ["0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],
  "2. Empty Line"     : ["\t\t\t\t","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],
  "3. Wrong Type"     : ["0.0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],
  "4. Invalid 'itra'" : ["-1 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],
  "5. Invalid 'iver'" : ["0 0.0 0.0 0.0 -1","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],
  "6. Too Many Lines" : ["0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9","Oops"],
}

with open("fort.3.template","r") as inFile:
  tmpF3 = inFile.read()

outBuf  = "%-32s %s\n" % ("Test Name","Exit=0")
outBuf += ("="*40)+"\n"
for aTest in sorted(theTests.keys()):
  theBlock = "\n".join(theTests[aTest][0])
  with open("fort.3","w") as outFile:
    outFile.write(tmpF3.replace("%INIT%",theBlock))
  exCode  = os.system(os.path.join("..","..","sixtrack"))
  outBuf += "%-32s %s\n" % (aTest,str(exCode == 0))

with open("error_results.log","w") as outFile:
  outFile.write(outBuf)
