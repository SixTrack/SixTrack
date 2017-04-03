#!/usr/bin/env python

import sys
import os,shutil
import re

CUTSPACE=True #Cut trailing spaces; old astuce seems to
              # - remove trailing spaces except for
              # - Add a single space to empty lines

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

def doWrite(ofile,line,level,tag=""):
        if inSUB or inFUN or inPRO:
            if CUTSPACE:
                line = line.rstrip()+'\n'
            ofile.write(str(level) + ":" + line)
        else:
            print "JUNK: "+line[:-1]+tag
            junklines.append(line)

isCONT = False
blocname = None
ofile = None
level=0
atSubEnd=False
inIF_cont = False
FUN_counter = 0
for line in ifile.xreadlines():
    if line[0] != ' ' and not line[0].isdigit():
        #Comment line -> skip
        doWrite(ofile,line,level)
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
            doWrite(ofile,line,level)
            continue
    print line[:-1], line_stripped
        
    #Look for the start of a new block
    if line_stripped.startswith("subroutine") or line_stripped.startswith("SUBROUTINE"):
        #print line[:-1]
        assert(inSUB==False and inFUN==False and inPRO==False), line
        inSUB = True
        
        assert blocname == None
        blocname = line_stripped.split()[1].split("(")[0].strip()

        print "New subrotine block:", blocname
        ofile = open(os.path.join(DIRNAME,"SUB-"+blocname+".f"),'w')

        level += 1
    elif (line_stripped.startswith("do") or line_stripped.startswith("DO"))\
         and not (line_stripped.startswith("double") or line_stripped.startswith("DOUBLE")):
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
    elif line_stripped.startswith("if") or line_stripped.startswith("IF"):
        #TODO: if statements may be short form, i.e. no "then".
        #print "if:",line_stripped
        if ("then" in line_stripped) or ("THEN" in line_stripped):
            level +=1
        else:
            #print "setting inIF_cont"
            inIF_cont = True # The "then" may come later...
    elif line_stripped.startswith("select") or line_stripped.startswith("SELECT"):
        level += 1
    elif (line_stripped.startswith("end") or line_stripped.startswith("END"))\
         and (not (line_stripped.startswith("endfile") or line_stripped.startswith("ENDFILE"))):
        #print "END: ", line[:-1]
        level -=1
        if level == 0:
            atSubEnd = True
    elif ("function" in line_stripped) or ("FUNCTION" in line_stripped):
        if ("'" in line_stripped) or ('"' in line_stripped):
            doWrite(ofile,line,level,"(FAKEFUN)")
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

    doWrite(ofile,line,level,"(OUTSIDE)")
        
    if atSubEnd:
        print "End of block",blocname
        ofile.close()
        blocname=None
        if inSUB:
            inSUB = False
        elif inFUN:
            inFUN = False
        atSubEnd = False
