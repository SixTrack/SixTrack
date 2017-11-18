#!/usr/bin/env python

import sys

assert len(sys.argv) == 3, "Usage: to_replace filename"

to_replace = sys.argv[1]
filename = sys.argv[2]

lines_in = open(filename,'r').readlines()

#print '+ca '+to_replace

i = 0
num_replaced = 0
while True:
    line = lines_in[i]
    #print line
    
    if line.startswith('+ca '+to_replace):
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
                break
        lines_in.insert(i,"      use "+to_replace + "\n")
        num_replaced += 1
    i = i+1
    if i >= len(lines_in):
        break

file_out = open(filename,'w')
for l in lines_in:
    file_out.write(l)
file_out.close()

print num_replaced
