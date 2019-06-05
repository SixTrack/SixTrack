import numpy as np
from ctypes import *
import matplotlib.pyplot as plt
from distpyinterface import *
from mpl_toolkits.mplot3d import Axes3D
def readtasfromdump(filename):
	file = open(filename, "r")
	splitlines = (file.readlines())
	tasstr = (str.split(splitlines[3]))
	tasstr = np.array(tasstr[3:])
	tas = tasstr.astype(np.float)
	tasm = tas.reshape(6,6)
	return tasm

def settasmatrix(dist, tas):
	for i in range(0,6):
		for j in range(0,6):
			dist.settasmatrix_element(c_double(tas[i][j]), c_int(i), c_int(j))


dist = cdll.LoadLibrary("./buildDemo/libhello.so")

dist.initializedistribution_(byref(c_int(2)))
myfile = "/home/tobias/codes/SixTrackTobias/test/orbit6d-element-quadrupole/START_DUMP"
tas = readtasfromdump(myfile)
settasmatrix(dist, tas)


pia2 = np.pi*2
momentum = c_double(6500)
mass = c_double(938.0)
eps = 2.0
dim = c_int(6)
zero = c_double(0)
e1 = 1.7
e2 = 1.0
e3 = 0.3

setEmittance12(dist,e1,e2)
setEmittance3(dist, e3)
setmassmom(dist, mass, momentum)
setdisttype(dist,1)

dist.settotallength(byref(c_int(5000)))
setParameters(dist,1,0,1,200,6)
setParameters(dist,2,0,pia2,200,4)
setParameters(dist,3,0,1,200,6)
setParameters(dist,4,0,pia2,200,4)
setParameters(dist,5,0,1,200,6)
setParameters(dist,6,0,pia2,200,4)


dist.printdistsettings_()
dist.print2file()


[xd,pxd,yd,pyd,deltas,dp] = getnormalizedcoordinates(dist)

x_n  = np.array(xd)
xp_n = np.array(pxd)
y_n  = np.array(yd)
yp_n = np.array(pyd)
d_n  = np.array(deltas)
dp_n = np.array(dp)


e1 = x_n**2+xp_n**2
e2 = y_n**2+yp_n**2
e3 = d_n**2+1000*dp_n**2

print(np.mean((e1)))
print(np.mean((e2)))
print(np.mean((e3)))

[xd,pxd,yd,pyd,deltas,dp] = getphysicalocordinates(dist)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(xd, yd, deltas)

plt.figure(2)
plt.scatter(xd,pxd)
#plt.hist(e3)

plt.show()