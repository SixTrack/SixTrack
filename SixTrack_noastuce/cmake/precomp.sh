#!/bin/bash

# input:
# $1: compiler
# $2: include flags
# $3: input file
# $4: output file
$1 $2 -E $3 | grep -v "^#" > $4
