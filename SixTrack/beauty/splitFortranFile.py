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
inIF_cont = False
FUN_counter = 0
for line in ifile.xreadlines():
    if line[0] != ' ':
        #Comment line -> skip
        if inSUB or inFUN or inPRO:
            ofile.write(str(level) + ":" + line)
        else:
            print "JUNK: "+line[:-1]
            junklines.append(line)
        continue
    
    #Don't consider labels!
    line_stripped = line[6:]
    #Don't consider anything past the end of the CARD
    CARDLEN = 66
    if len(line_stripped) > CARDLEN: 
        line_stripped=line_stripped[:CARDLEN]
    #Don't consider the part of the line that is a comment
    commentIdx = line_stripped.find('!')
    if commentIdx != -1:
        line_stripped = line_stripped[:commentIdx]
    # #Don't consider whatever is inside a string
    # startStrStrip=0
    # while True:
    #     strIdx1 = line_stripped.find("'",startStrStrip)
    #     if strIdx1 != -1:
    #         strIdx2 = line_stripped.find("'",strIdx1+1)
    #         print line_stripped
    #         if strIdx2 == -1: #String runs to the end of the line
    #             print "type1:"
    #             print line_stripped[:-1]
    #             line_stripped=line_stripped[:strIdx1+1]
    #             print line_stripped[:-1]
    #             startStrStrip = -1
    #         else: #String runs to the end of the line
    #             print "type2:"
    #             print line_stripped[:-1]
    #             line_stripped=line_stripped[:strIdx1+1]+line_stripped[strIdx2:]
    #             print line_stripped[:-1]
    #             startStrStrip = strIdx2+1
    #     else:
    #         #print "Done."
    #         break
    
    #Identify continuation lines
    isCONT = False
    if len(line)>5 and line[5] != ' ':
        isCONT = True
    else:
        inIF_cont = False
    line_stripped = line_stripped.strip()
    #Some if ... then with continuation
    if isCONT and inIF_cont:
        if "then" in line_stripped:
            level +=1
            if inSUB or inFUN or inPRO:
                ofile.write(str(level) + ":" + line)
            else:
                print "JUNK: "+line[:-1]
                junklines.append(line)
            continue

    #Look for the start of a new block
    if line_stripped.startswith("subroutine"):
        #print line[:-1]
        assert(inSUB==False and inFUN==False and inPRO==False), line
        inSUB = True
        
        assert blocname == None
        blocname = line_stripped.split()[1].split("(")[0].strip()

        print "New subrotine block:", blocname
        ofile = open(os.path.join(DIRNAME,"SUB-"+blocname+".f"),'w')

        level += 1
    elif line_stripped.startswith("do") and not line_stripped.startswith("double"):
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
            #print "setting inIF_cont"
            inIF_cont = True # The "then" may come later...
    elif line_stripped.startswith("select"):
        level += 1
    elif line_stripped.startswith("end") and (not line_stripped.startswith("endfile")):
        #print "END: ", line[:-1]
        level -=1
        if level == 0:
            atSubEnd = True
    elif "function" in line_stripped:
        if ("'" in line_stripped) or ('"' in line_stripped):
            if inSUB or inFUN or inPRO:
                ofile.write(str(level) + ":" + line)
            else:
                print "JUNK: "+line[:-1]+ "(FAKEFUN)"
                junklines.append(line)
        else:
            #print line[:-1]
            assert(inSUB==False and inFUN==False and inPRO==False), (inSUB, inFUN, inPRO, line)
            inFUN = True
        
            assert blocname == None
            #blocname = line_stripped.split()[1].split("(")[0].strip()
            blocname=str(FUN_counter)
            FUN_counter += 1
            
            print "New function block:", blocname
            ofile = open(os.path.join(DIRNAME,"FUN-"+blocname+".f"),'w')

            level += 1
    assert level >= 0

    if inSUB or inFUN or inPRO:
        ofile.write(str(level) + ":" + line)
    else:
        print "JUNK: "+line[:-1] + "(OUTSIDE)"
        print line_stripped[:-1]
        junklines.append(line)
        #exit (1)
        
    if atSubEnd:
        print "End of block",blocname
        ofile.close()
        blocname=None
        if inSUB:
            inSUB = False
        elif inFUN:
            inFUN = False
        atSubEnd = False
