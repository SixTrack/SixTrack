#!/usr/bin/env python

import sys, os

defLn1 = (0,0.0,0.0,0.0,0)
defLnN = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,7.0e9,7.0e9,7.0e9)

theTests = {
  "Working INIT"       : [["0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],None],
  "Ln1 Empty"          : [["\t\t\t\t","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],None],
  "Ln1 Wrong Type"     : [["0.0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],None],
  "Ln1 Invalid 'itra'" : [["-1 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],None],
  "Ln1 Invalid 'iver'" : [["0 0.0 0.0 0.0 -1","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9"],None],
  "Too Many Lines"     : [["0 0.0 0.0 0.0 0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","7.0e9","7.0e9","7.0e9","Oops"],None],
}

with open("fort.3.template","r") as inFile:
  tmpF3 = inFile.read()

for aTest in theTests.keys():
  theBlock = "\n".join(theTests[aTest][0])
  with open("fort.3","w") as outFile:
    outFile.write(tmpF3.replace("%INIT%",theBlock))
  theTests[aTest][1] = os.system(os.path.join("..","..","sixtrack"))

with open("error_results.log","w") as outFile:
  outFile.write("%-32s %s\n" % ("Test Name","Exit=0"))
  outFile.write(("="*40)+"\n")
  for aTest in theTests.keys():
    if theTests[aTest][1] is None:
      exStat = None
    else:
      exStat = theTests[aTest][1] == 0
    outFile.write("%-32s %s\n" % (aTest,str(exStat)))
