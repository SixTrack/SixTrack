#!/usr/bin/env python
import numpy as np

#npart = 1936 #BIGNPART
#npart = 256**2 # HUGENPART, use this to make nice plots
npart = 64**2 # HUGENPART, reduced for the sake of SixTest
assert npart % 2 == 0
npart_1d = int(np.sqrt(npart))
assert npart_1d**2 == npart

## BEAM PARAMETERS

#Rectangular grid, covering ~5 sigma of the bb beam
xmax= 3*np.sqrt( 2.1046670129999999e-05 ) #[mm]
ymax= 3*np.sqrt( 3.165637487e-06        ) #[mm]

h_x = 2*xmax/float(npart_1d-1)
h_y = 2*ymax/float(npart_1d-1)
print "hx,hy=",h_x,h_y,"[mm]"

E0    = 0.45e12   # Beam energy          [eV]
epsx_n = 2.5e-6   # Normalized emittance [m*rad]
epsy_n = 2.5e-6   # Normalized emittance [m*rad]

#Closed orbit
x0 = 0.0
y0 = 0.0
xp0 = 0
yp0 = 0.0

#General physics parameters
mp = 938.272046e6 #proton mass, eV/c^2
e  = 1.60217657e-19 #C, electron charge
c  = 2.99792485e8   #m/s, speed of light

#Calculate parameters
gamma_rel = E0/mp
beta_rel  = np.sqrt(1-gamma_rel**-2)
p0 = np.sqrt((E0-mp)*(E0+mp))

x   = np.zeros(npart)
xp  = np.zeros(npart)
y   = np.zeros(npart)
yp  = np.zeros(npart)
z   = np.zeros(npart)
dPP = np.zeros(npart)
E   = np.zeros(npart)

i = 0
for i_x in xrange(npart_1d):
    for i_y in xrange(npart_1d):
        x[i]  = -xmax + i_x*h_x + x0
        xp[i] = 0.0 + xp0
        y[i]  = -ymax + i_y*h_y + y0
        xp[i] = 0.0 + yp0

        z[i]   = 0.0
        dPP[i] = 0.0
        E[i]   = E0
        
        i+=1

#Write to file
ofile = open("fort.13", 'w')
for i in xrange(0,npart,2):
    ofile.write(str(x  [i] )+"\n") #mm
    ofile.write(str(xp [i] )+"\n") #mrad
    ofile.write(str(y  [i] )+"\n") #mm
    ofile.write(str(yp [i] )+"\n") #mrad
    ofile.write(str(z  [i] )+"\n") #mm
    ofile.write(str(dPP[i] )+"\n") #-

    ofile.write(str(x  [i+1] )+"\n") #mm
    ofile.write(str(xp [i+1] )+"\n") #mrad
    ofile.write(str(y  [i+1] )+"\n") #mm
    ofile.write(str(yp [i+1] )+"\n") #mrad
    ofile.write(str(z  [i+1] )+"\n") #mm
    ofile.write(str(dPP[i+1] )+"\n") #-

    ofile.write(str(E0*1e-6      )+"\n") #MeV
    ofile.write(str(E  [i]  *1e-6)+"\n") #MeV
    ofile.write(str(E  [i+1]*1e-6)+"\n") #MeV
ofile.close()

import matplotlib.pyplot as plt
plt.plot(x,y,'+')
plt.show()
