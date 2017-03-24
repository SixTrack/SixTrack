#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print "Usage: ./"+sys.argv[0]+" FILENAME.s"
    
ifile = open(sys.argv[1],'r')
ofile = open("COUNTED-"+sys.argv[1],'w')
ifcounter = 0
linecounter = 0
for l in ifile.xreadlines():
    change = False
    linecounter += 1
    if l.startswith("+if"):
        ifcounter += 1
        change = True
    elif l.startswith("+ei"):
        ifcounter -= 1
        change = True

#    if change==True and ifcounter <= 0:
#        print "ifcounter=",ifcounter," at line=",linecounter

    ofile.write(str(ifcounter) + " " + l)
    
ifile.close()
ofile.close()
print "Done! Please see file 'COUNTED-"+sys.argv[1]+"' !"
