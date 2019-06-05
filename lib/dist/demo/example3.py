from ctypes import *
import numpy
import matplotlib.pyplot as plt
from distpyinterface import *


def setParameters(dist, index, start, stop, numb, type):
	dist.setparameter_(byref(c_int(index)), byref(c_double(start)),byref(c_double(stop)),byref(c_int(numb)),byref(c_int(type)))


double6 = c_double * 6
DOUBLE = c_double
PDOUBLE = POINTER(DOUBLE)

eps = 2.0
dim = c_int(6)
zero = c_double(0)
e1 = (eps)
e2 = (1.0)
e3 = (1.0)
dp = c_double(0.00027)
pia2 = numpy.pi*2
zero = 0
one = 1
num = 10000
betx=1
alfx= -0.50
bety=1
alfy= 0

momentum = c_double(6500000)
mass = c_double(938.0)
dim = c_int(6)
coord_c = double6(1,1,1,1,1,0)
dist = cdll.LoadLibrary("./buildDemo/libhello.so")

thdeg = numpy.linspace(0,2*numpy.pi, 100)
teta =[]
xa = []
xpa = []
beta= betx
alfa =alfx
for i in range(0, len(thdeg)):
	teta.append(thdeg[i])
	xa.append(numpy.sqrt(eps*beta)*numpy.cos(teta[i]))
	xpa.append(-numpy.sqrt(eps/beta)*( alfa*numpy.cos(teta[i]) + numpy.sin(teta[i]) ))

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

acoord = double6(1,1,1,1,0,0)
physical = double6(0,0,0,0,0,0)


dist.initializedistribution_(byref(c_int(2)))
setEmittance12(dist,e1,e2)
setEmittance3(dist, e3)
setmassmom(dist, mass, momentum)
setdisttype(dist,0)
createtas0coupling(dist, betx,alfx,bety,alfy, zero, zero, zero, zero)


setParameters(dist,1,1,4,4,1)
setParameters(dist,2,0,pia2,200,1)
setParameters(dist,3,1,4,1,1)
setParameters(dist,4,0,0,1,0)
setParameters(dist,5,0.000,0,1,0)
setParameters(dist,6,0.000,0,1,0)

dist.setnormalizedcut(c_int(1),c_double(0), c_double(100))
dist.setphysicalcut(c_int(2),c_double(0), c_double(100))


[x,px,y,py,deltas,dp] = getphysicalocordinates(dist)

print(dist.getnumberdist_())

plt.plot(x,px, '*')
plt.plot(xa,xpa, '*')

dist.printdistsettings_()
plt.show()

#dim_c = c_int(6)
#dist.addclosedorbit_(coord_c)
#mychar = c_char_p(b"wrrritiing")
#dist.printvector(mychar, dim_c, coord_c)
#for d in coord_c:
#	print(d)
