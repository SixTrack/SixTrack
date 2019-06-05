from ctypes import *
import numpy
import matplotlib.pyplot as plt
from distpyinterface import *


eps=1
e1 = (eps)
e2 = (1.0)
e3 = (1.0)
#dp =_double(0.00027)
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


dist.initializedistribution_(byref(c_int(2)))
setEmittance12(dist,e1,e2)
setEmittance3(dist, e3)
setmassmom(dist, mass, momentum)
setdisttype(dist,0)
createtas0coupling(dist, betx,alfx,bety,alfy, zero, zero, zero, zero)

setParameters(dist,1,2,0,0,0	)
setParameters(dist,2,0,pia2,200,1)
setParameters(dist,3,1,4,1,1)
setParameters(dist,4,0,0,1,0)
setParameters(dist,5,0.001,0,1,0)
setParameters(dist,6,0,0,1,0)


[x,px,y,py,deltas,dp] = getphysicalocordinates(dist)

print(dist.getnumberdist_())

plt.plot(x,px, '*')
plt.plot(xa,xpa, '*')

dist.printdistsettings_()
plt.show()

