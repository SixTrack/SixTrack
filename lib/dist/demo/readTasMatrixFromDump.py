import numpy as np
from ctypes import *
import matplotlib.pyplot as plt
from distpyinterface import *

def readtasfromdump(filename):
	file = open(filename, "r")
	splitlines = (file.readlines())
	tasstr = (str.split(splitlines[3]))
	tasstr = np.array(tasstr[3:])
	tas = tasstr.astype(np.float)
	tasm = tas.reshape(6,6)
	return tasm

def setParameters(dist, index, start, stop, numb, type):
	dist.setparameter_(byref(c_int(index)), byref(c_double(start)),byref(c_double(stop)),byref(c_int(numb)),byref(c_int(type)))


dist = cdll.LoadLibrary("./buildDemo/libhello.so")


myfile = "/home/tobias/codes/SixTrackTobias/test/orbit6d-element-quadrupole/START_DUMP"
tas = readtasfromdump(myfile)
settasmatrix(dist, tas)



double6 = c_double * 6
physical = double6(0,0,0,0,0,0)
pia2 = np.pi*2
momentum = c_double(6500000)
mass = c_double(938.0)
eps = 2.0
dim = c_int(6)
zero = c_double(0)
e1 = c_double(1.0)
e2 = c_double(1.0)
e3 = c_double(0.0)
dist.initializedistribution_(byref(c_int(2)))
dist.setemittance12_(byref(e1),byref(e2))
dist.setemittance3_(byref(e3))
dist.setmassmom_(byref(mass), byref(momentum))
dist.setdisttype(byref(c_int(1)))	


dist.settotallength(byref(c_int(2000)))

setParameters(dist,1,0,1,200,6)
setParameters(dist,2,0,pia2,200,4)
setParameters(dist,3,0,1,200,6)
setParameters(dist,4,0,0,1,0)
setParameters(dist,5,0,0,1,0)
setParameters(dist,6,0,0,1,0)

dist.printdistsettings_()
dist.print2file()

xd = []
pxd = []
yd = []
yxd = []


for i in range(0,dist.getnumberdist_()):
	dist.getcoord_(physical,c_int(i))
	xd.append(physical[0])
	pxd.append(physical[1])
	yd.append(physical[2])
	yxd.append(physical[3])
	


dist.setdistribution_(byref(c_int(1)))
dist.setemittance12_(byref(e1),byref(e2))
dist.setemittance3_(byref(e3))
dist.setmassmom_(byref(mass), byref(momentum))
dist.setdisttype(byref(c_int(0)))
settasmatrix(dist, tas)


dist.setdisttype(byref(c_int(0)))


setParameters(dist,1,1,1,1,1)
setParameters(dist,2,0,pia2,200,1)
setParameters(dist,3,1,1,1,1)
setParameters(dist,4,0,0,1,0)
setParameters(dist,5,0.000,0,1,0)
setParameters(dist,6,0,0,1,0)

xd_c = []
pxd_c = []
yd_c = []
yxd_c = []

for i in range(0,dist.getnumberdist_()):
	dist.getcoord_(physical,c_int(i))
	xd_c.append(physical[0])
	pxd_c.append(physical[1])
	yd_c.append(physical[2])
	yxd_c.append(physical[3])
	


print(tas)
print(len(xd_c))
dist.printdistsettings_()
dist.print2fort13()

plt.plot(xd, pxd, '*')
plt.plot(xd_c, pxd_c, '*')
plt.show()