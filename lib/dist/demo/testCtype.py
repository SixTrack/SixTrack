from ctypes import *
import numpy
import matplotlib.pyplot as plt
a = numpy.loadtxt('data/tasm.txt')
mat = a.reshape(6,6)


double6 = c_double * 6
DOUBLE = c_double
PDOUBLE = POINTER(DOUBLE)
ptr = (DOUBLE*6*6)()
DBL5ARR = DOUBLE * 6
# An array of double* can be passed to your function as double**.
PDBL4ARR = PDOUBLE * 6

ptr = PDBL4ARR()
for i in range(6):
    # fill out each pointer with an array of doubles.
    ptr[i] = DBL5ARR()
    for j in range(6):
        ptr[i][j] = mat[i][j]  # just to initialize the actual doubles.


dim = c_int(6)

coord_c = double6(1,1,1,1,1,0)
dist = cdll.LoadLibrary("./buildDemo/libhello.so")
#dist.hello()

#dist.printmatrix(dim, dim, ptr)

dim = c_int(6)
e1 = c_double(1.0)
e2 = c_double(1.0)
e3 = c_double(1.0)
dp = c_double(0.0000)
pia2 = c_double(numpy.pi*2)
betx=1
alfx= 0
bety=1
alfy= 0
momentum = c_double(4000)
mass = c_double(938.0)


acoord = double6(1,1,1,1,0,0)
physical = double6(0,0,0,0,0,0)

dist.initializedistribution_(byref(dim), byref(dim))
print("initial")
# Set the tas matrix 
#dist.settasmatrixpython(ptr)
dist.createtas0coupling_(c_double(betx),c_double(alfx),c_double(bety),c_double(alfy))
# Set the emittance
dist.setemittance12_(byref(e1),byref(e2))
dist.setemittance3_(byref(dp))
#dist.setemittance3(e3)
dist.setmassmom_(byref(mass), byref(momentum))
dist.action2canonical_(acoord,physical)



thdeg = numpy.linspace(0,2*numpy.pi, 100)
teta =[]
x = []
xp = []
beta= betx
alfa =alfx
eps = 3*1#2*numpy.pi
for i in range(0, len(thdeg)):
	teta.append(thdeg[i])
	x.append(numpy.sqrt(eps*beta)*numpy.cos(teta[i]))
	xp.append(-numpy.sqrt(eps/beta)*( alfa*numpy.cos(teta[i]) + numpy.sin(teta[i]) ))




xa = []
y = []
for i in range(0,10000):
	dist.createrandom(acoord, physical)
	xa.append(physical[0])
	y.append(physical[1])

#plt.hist(x,bins=50)
plt.plot(xa,y, '.')
plt.plot(x,xp)
print(numpy.mean(numpy.array(xa)**2))
print(numpy.mean(numpy.array(xa)*numpy.array(y)))
print(numpy.std(numpy.array(xa)))
count = 0 
for i in range(0,10000):
	if(numpy.sqrt(xa[i]**2+y[i]**2) <3):
		count = count+1
print(count/10000)
plt.show()
#plt.hist(y,bins=50)
#plt.show()


#dim_c = c_int(6)
#dist.addclosedorbit_(coord_c)
#mychar = c_char_p(b"wrrritiing")
#dist.printvector(mychar, dim_c, coord_c)
#for d in coord_c:
#	print(d)
