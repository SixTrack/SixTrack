#!/usr/bin/env python

import sys
import os,shutil
import re

if len(sys.argv) != 2:
    print "Usage:"
    print sys.argv[0] + " file_to_split.f"
    exit(1)

ifile = open(sys.argv[1], 'r')

DIRNAME = sys.argv[1]+"-SPLITTED"

if os.path.exists(DIRNAME):
    if not os.path.isdir(DIRNAME):
        print "ERROR: '"+DIRNAME+"' exists but is not a directory."
        exit(1)
    shutil.rmtree(DIRNAME)
os.mkdir(DIRNAME)

junklines = []

inSUB = False
inFUN = False
inPRO = False
isCONT = False
blocname = None
ofile = None
level=0
atSubEnd=False
for line in ifile.xreadlines():
    isCONT = False
    #Don't consider labels!
    line_stripped = line[6:]
    #Don't consider anything past the end of the CARD
    CARDLEN = 66
    if len(line_stripped) > CARDLEN: 
        line_stripped=line_stripped[:CARDLEN]
    if len(line)>5 and line[5] != ' ':
        isCONT = True
    else:
        inIF_cont = False
    line_stripped = line_stripped.strip()
    if line[0] != ' ':
        #Comment line -> skip
        continue
    
    if isCONT and inIF_cont:
        if "then" in line_stripped:
            level +=1
    
    if line_stripped.startswith("subroutine"):
        #print line[:-1]
        assert(inSUB==False and inFUN==False and inPRO==False), line
        inSUB = True
        
        assert blocname == None
        blocname = line_stripped.split()[1].split("(")[0].strip()

        print "New subrotine block:", blocname
        ofile = open(os.path.join(DIRNAME,"SUB-"+blocname+".f"),'w')

        level += 1
    elif line_stripped.startswith("do"):
        do_tmp = line_stripped[2:].strip()
        if do_tmp[0].isdigit():
            #old style 'do LABEL'
            #print "dolabel:", line_stripped
            pass
        elif do_tmp[:4]=="uble":
            #double precision something something...
            #print "double:", line_stripped
            pass
        else:
            #print line[:-1]
            #print "do:",line_stripped
            level +=1
    elif line_stripped.startswith("if"):
        #TODO: if statements may be short form, i.e. no "then".
        #print "if:",line_stripped
        if "then" in line_stripped:
            level +=1
        else:
            inIF_cont = True # The "then" may come later...
    elif line_stripped.startswith("select"):
        level += 1
    elif line_stripped.startswith("end"):
        #print line[:-1]
        level -=1
        if level == 0:
            atSubEnd = True
        
    assert level >= 0

    if inSUB or inFUN or inPRO:
        ofile.write(str(level) + " " + line)
    else:
        print "JUNK: "+line[:-1]
        junklines.append(line)
        #exit (1)
        
    if atSubEnd:
        print "End of block",blocname
        ofile.close()
        blocname=None
        if inSUB:
            inSUB = False
        atSubEnd = False
