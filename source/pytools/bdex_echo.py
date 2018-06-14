#!/usr/bin/env python

# ECHO particles for BDEX
# This reads the particles from BDEX, and writes back the same particles.

import sys

if len(sys.argv) != 3:
    print "Arguments: InPipe OutPipe"
    exit(1)

oPipeName = sys.argv[1] #/tmp/pip3
iPipeName = sys.argv[2] #/tmp/pip4

#COPYBACK=False
COPYBACK=True

print "Opening", oPipeName, "for writing"
oPipe = open(oPipeName,'w')
print "Opening", iPipeName, "for reading"
iPipe = open(iPipeName,'r')

isConnected=False

while True:
    iData = iPipe.readline()
    if len(iData)==0:
        print "No more input!"
        break
    iData = iData.strip()
    
    if iData.startswith("BDEX-PIPE"):
        print "SixTrack connected."
        assert isConnected==False
        isConnected=True
    elif iData.startswith("BDEX TURN"):
        ls = iData.split()
        turn = int(ls[2])
        bez  = ls[3][4:]
        I    = int(ls[5])
        NAPX = int(ls[7])
        particles = []
        for j in xrange(NAPX):
            particles.append( iPipe.readline().strip() )
            #print particles[-1]
        
        waiting=iPipe.readline().strip()
        #print waiting
        assert waiting=="BDEX WAITING..."
        print "Got", NAPX, "particles"
        
        if (COPYBACK==True):
            oPipe.write(str(NAPX)+"\n")
            for j in xrange(NAPX):
                oPipe.write(particles[j]+"\n")
            oPipe.flush()
            print "Wrote them back!"
        else:
            oPipe.write("-1\n")
            oPipe.flush()
            print "Wrote -1 -> No change"
        
        tracking=iPipe.readline().strip()
        #print tracking
        assert tracking=="BDEX TRACKING..."
        
iPipe.close()
oPipe.close()
