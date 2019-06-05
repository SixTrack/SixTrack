from ctypes import *
def setParameters(dist, index, start, stop, numb, type):
	dist.setparameter_(byref(c_int(index)), byref(c_double(start)),byref(c_double(stop)),byref(c_int(numb)),byref(c_int(type)))

def setEmittance12(dist, e1, e2):
	dist.setemittance12_(byref(e1),byref(e2))

def setEmittance3(dist, e3):
	dist.setemittance3_(byref(e3))

def createtas0coupling(dist, betx, alfx, bety, alfy, dx, dpx, dy, dpy):
	dist.createtas0coupling_(c_double(betx),c_double(alfx),c_double(bety),c_double(alfy), c_double(0), c_double(0), c_double(0), c_double(0))

def setmassmom(dist, mass, momentum):
	dist.setmassmom_(byref(mass), byref(momentum))

def setdisttype(dist, disttype):
	dist.setdisttype(byref(c_int(disttype)))

def getphysicalocordinates(dist):
	double6 = c_double * 6
	physical = double6(0,0,0,0,0,0)
	xd = []
	pxd = []
	yd = []
	yxd = []
	deltas = []
	dp = []
	for i in range(0,dist.getnumberdist_()):
		dist.getcoord_(physical,c_int(i))
		xd.append(physical[0])
		pxd.append(physical[1])
		yd.append(physical[2])
		yxd.append(physical[3])
		deltas.append(physical[4])
		dp.append(physical[5])
	return [xd, pxd, yd, yxd, deltas, dp]