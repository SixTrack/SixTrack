#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

quiverSkip=8

# Read the dump file from just after the bb kick

fileDType = np.dtype([('ID', np.int), ('turn', np.int),('s', np.float),
                      ('x', np.float),('xp', np.float),
                      ('y', np.float),('yp', np.float),
                      ('z', np.float),('dEE', np.float),
                      ('ktrack', np.int)])
dumpfile = np.loadtxt("dump_bb.dat",dtype=fileDType)
dumpfile_initial = np.loadtxt("dump_ip.dat",dtype=fileDType)

#
npart = len(dumpfile['x'])
npart_1d = int(np.sqrt(npart))
assert npart_1d**2 == npart

x = dumpfile['x']
x.shape=(npart_1d,npart_1d)
y = dumpfile['y']
y.shape=(npart_1d,npart_1d)

xp = dumpfile['xp']
xp.shape=(npart_1d,npart_1d)
yp=dumpfile['yp']
yp.shape=(npart_1d,npart_1d)

x_initial = dumpfile_initial['x']
x_initial.shape=(npart_1d,npart_1d)
y_initial = dumpfile_initial['y']
y_initial.shape=(npart_1d,npart_1d)

h_x = x_initial[1,0]-x_initial[0,0]
h_y = y_initial[0,1]-y_initial[0,0]
print "hx,hy=",h_x,h_y,"[mm]"

#Calculate rho, proportional to the charge density of the strong beam
DxpDx = np.zeros((npart_1d-2,npart_1d-2))
DypDy = np.zeros((npart_1d-2,npart_1d-2))
for ix in xrange(1,npart_1d-1):
    for iy in xrange(1,npart_1d-1):
        DxpDx[ix-1,iy-1] = (xp[ix+1,iy]-xp[ix-1,iy])/(2*h_x)
        DypDy[ix-1,iy-1] = (yp[ix,iy+1]-yp[ix,iy-1])/(2*h_y)
rho = DxpDx + DypDy

maxKick=np.max(np.sqrt(xp**2+yp**2))
print "maxKick=",maxKick,"mrad"

#Kicks and charge density
plt.figure()

if maxKick>0.0:
    plt.contour(x_initial[1:-1,1:-1],y_initial[1:-1,1:-1],rho,20)
    plt.colorbar()

plt.quiver(x_initial[::quiverSkip,::quiverSkip],y_initial[::quiverSkip,::quiverSkip],xp[::quiverSkip,::quiverSkip],yp[::quiverSkip,::quiverSkip])

plt.xlabel("x [mm]")
plt.ylabel("y [mm]")

#Plot as a function of x,y on contours running through (0,0)
#Find index of x=0,y=0:
xi_0 = -1
x_min = 1000.0
yi_0 = -1
y_min = 1000.0

assert x_initial.shape[0]==x_initial.shape[1]
for i in xrange(1,x_initial.shape[0]-1):
    if abs(x[i+1,0]) < x_min:
        x_min = abs(x[i+1,0])
        xi_0 = i
    if abs(y[0,i+1]) < y_min:
        y_min = abs(y[0,i+1])
        yi_0 = i
print x_min,xi_0
print y_min,yi_0

plt.axhline(y_min,ls="--",color="r")
plt.axvline(x_min,ls="--",color="r")

plt.figure()
plt.plot(x[1:-1,yi_0+1],rho[:,yi_0],label="Hor.")
plt.plot(y[xi_0+1,1:-1],rho[xi_0,:],label="Ver.")
plt.legend()
plt.xlabel("Distance along contour [mm]")

#Plot Delta x, Delta y:
plt.figure()

xShift = x-x_initial
yShift = y-y_initial
maxShift=np.max(np.sqrt(xShift**2+yShift**2))

plt.quiver(x_initial[::quiverSkip,::quiverSkip],y_initial[::quiverSkip,::quiverSkip],xShift[::quiverSkip,::quiverSkip],yShift[::quiverSkip,::quiverSkip])
plt.xlabel("x [mm]")
plt.ylabel("y [mm]")
plt.title("Point displacement; max = %g [mm]" %(maxShift,))

plt.show()
