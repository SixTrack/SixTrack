#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# Read the dump file from just after the bb kick

fileDType = np.dtype([('ID', np.int), ('turn', np.int),('s', np.float),
                      ('x', np.float),('xp', np.float),
                      ('y', np.float),('yp', np.float),
                      ('z', np.float),('dEE', np.float),
                      ('ktrack', np.int)])
dumpfile = np.loadtxt("dump_bb.dat",dtype=fileDType)

plt.quiver(dumpfile["x"],dumpfile["y"],dumpfile["xp"],dumpfile["yp"])
plt.show()
