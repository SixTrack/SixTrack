#!/usr/bin/env python

import sys

assert len(sys.argv) == 3 or len(sys.argv) == 4, "Usage: to_replace filename (modname)"

to_replace = sys.argv[1]
filename = sys.argv[2]
if len(sys.argv) == 4:
    modname = sys.argv[3]
else:
    modname = to_replace

string_replace = '+ca '+to_replace

lines_in = open(filename,'r').readlines()

#print '+ca '+to_replace

i = 0
num_replaced = 0
numspaces = 6 # default value
while True:
    line = lines_in[i]
    #print line
    
    if line.startswith(string_replace) and (line[len(string_replace)] == "\n" or line[len(string_replace)] == " "):
        #print "found!"
        #delete the bad line
        del lines_in[i]
        
        #search backwards for the implicit none
        while True:
            i = i-1
            if i < 0:
                print "Error, i<0"
                exit(1)
            line = lines_in[i]
            if "implicit none" in line or "IMPLICIT NONE" in line:
                numspaces = len(line)-len(line.lstrip())
                break
        lines_in.insert(i,numspaces*" " + "use " + modname + "\n")
        num_replaced += 1
    i = i+1
    if i >= len(lines_in):
        break

file_out = open(filename,'w')
for l in lines_in:
    file_out.write(l)
file_out.close()

print "num_replaced=", num_replaced
