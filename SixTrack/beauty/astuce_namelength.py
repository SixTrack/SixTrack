#!/usr/bin/env python

import sys

maxlen=8

if len(sys.argv) != 2:
    print "Checks for conflicts between flag/deck/block names"
    print "ASTUCE only cares about the first",maxlen,"characters."
    print "Usage: "+sys.argv[1]+" filename.s"
    exit(1)

    
flags = []
decks = []
blocks = []

allGood = True

ifile = open(sys.argv[1],'r')
lcounter = 0
for line in ifile.xreadlines():
    lcounter += 1
    if not line[0]=="+":
        continue
    elif line[:3] == "+ei":
        continue
    elif line[:3] == "+ca":
        continue
    
    if line[:3] == "+cd":
        ls = line.split()
        bname = ls[1]
        if bname in blocks:
            print "WARNING: Duplicate definition of block '"+bname+"'"
            print "Line: '"+line[:-1]+"'"
            print "Line number =", lcounter
            allGood=False
        else:
            blocks.append(bname)
    if line[:3] == "+if":
        flagsArray = []
        ls = line.split()
        ls1 = ls[1].split(".and.")
        for l in ls1:
            ls2=l.split(".or.")
            for ll in ls2:
                if ll[:5]==".not.":
                    flagsArray.append(ll[5:])
                else:
                    flagsArray.append(ll)
        #print line, flagsArray
        for flag in flagsArray:
            if len(flag) > 8:
                print "WARNING: Flag '"+ flag + "' too long!"
                allGood=False

    if line[:3] == "+dk":
        #print line
        ls = line.split()
        dname = ls[1]
        if len(dname) > 8:
            print "WARNING: Deck name '"+ dname + "' too long!"
            allGood=False

#print blocks
for i in xrange(len(blocks)):
    bi = blocks[i][:8]
    for j in xrange(i+1,len(blocks)):
        bj = blocks[j][:8]
        if bi==bj:
            print "WARNING: Blocks with duplicate first",maxlen,"characters",\
                "'"+bi+"': '"+blocks[i]+"' and '"+blocks[j]+"'"
            allGood=False

if not allGood:
    exit(1)
