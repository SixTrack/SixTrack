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

def settasmatrix(dist, tas):
	for i in range(0,6):
		for j in range(0,6):
			dist.settasmatrix_element(c_double(tas[i][j]), c_int(i), c_int(j))

def setParameters(dist, index, start, stop, numb, type):
	dist.setparameter_(byref(c_int(index)), byref(c_double(start)),byref(c_double(stop)),byref(c_int(numb)),byref(c_int(type)))

'''
SIMULATION
#/                   mass, charge,   A,    Z, momentum (7Z TeV/c)
REFERENCE 938.0      1  1    1  6500000.0 
/      turns, total particles
TRACK  100000   64
NEXT
/
DISTRIBUTION
//INPUT PARAMETERS
EMITTANCE   1           2.5e-6    //Normalized [m]
EMITTANCE   2           2.503e-6  //Normalized [m]
/          BETAX           ALFX            DX            DPX       BETAY          ALFY            DY              DPY
TWISS      156.921221293   2.11599999399   0.2973032402 -0.0030478 78.0330611129 -1.08169977617  -0.1155854323 -0.0011203
CLOSEDORBIT  0.0   0.0   0.0   0.0  0.0   0.0
DUMP FILENAME BINARY NORMALIZED
// DISTRIBUTION TYPE
SUBSAMPLE  1
JX       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
JY       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
DELTA    CONSTANT  0.00027
SIGMA    CONSTANT  0
SUBSAMPLE  2 /second pair
JX       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0.001
JY       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
DELTA    CONSTANT  0.00027
SIGMA    CONSTANT  0
NEXT
DUMP
POSTPR 1000 / produce fort.90 every 1000 turns
NEXT
'''

dist = cdll.LoadLibrary("./buildDemo/libhello.so")
dist.initializedistribution_(byref(c_int(2)))

myfile = "/home/tobias/codes/SixTrackTobias/test/orbit6d-element-dispersion/START_DUMP"
tas = readtasfromdump(myfile)
tas[4][2]=0.002
tas[5][2]=0.1
settasmatrix(dist, tas)



double6 = c_double * 6
physical = double6(0,0,0,0,0,0)
pia2 = np.pi*2
momentum = c_double(6500000)
mass = c_double(938.0)
eps = 2.0
dim = c_int(6)
zero = c_double(0)
e1 = (1.0)
e2 = (1.0)

setEmittance12(dist,e1,e2)
dist.usedeltap_()
setmassmom(dist, mass, momentum)
setdisttype(dist,0)


setParameters(dist,1,2,4,4,1)
setParameters(dist,2,0.001,0,1,0)
setParameters(dist,3,2,4,4,1)
setParameters(dist,4,0,0,1,0)
setParameters(dist,5,0.2,2,1,0)
setParameters(dist,6,1,1,1,0)

dist.printdistsettings_()
dist.print2file()

xd = []
pxd = []
yd = []
yxd = []
pyd = []
delta = []
time = []
for i in range(0,dist.getnumberdist_()):
	dist.getcoord_(physical,c_int(i))
	xd.append(physical[0])
	pxd.append(physical[1])
	yd.append(physical[2])
	pyd.append(physical[3])
	time.append(physical[4])
	delta.append(physical[5])

	print(physical[0], physical[1],physical[2],physical[3],physical[4], physical[5])


#dist.printdistsettings_()
#dist.print2fort13()

