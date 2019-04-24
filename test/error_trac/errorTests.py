#!/usr/bin/env python3

import sys
from errorTestTools import *

exCode = runTests({
  "0. Valid TRAC"     : ["1 0 1 0.0 0.0","1 1 1 1 1","0 0 1 1 1 50000 2"],
  "1. Empty Line"     : ["\t\t\t\t","1 0 1 0.0 0.0","1 1 1 1 1","0 0 1 1 1 50000 2"],
},sys.argv)

sys.exit(exCode)
