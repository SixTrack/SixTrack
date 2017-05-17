#!/usr/bin/env python
import numpy as np

import sys
# if len(sys.argv) ==2:
#     seed = int(sys.argv[1])
#     print "Got seed =", seed
#     np.random.seed(seed)
# else:
#     seed = 0

npart_amp = 8
npart_ang = 8
npart = npart_amp*npart_ang

amp_max = 8.1
amp_min = 0.1
ang_off = 0.0

#npart = 64   #Number of particles to generate
#             #Should be divisible by 2!

## BEAM PARAMETERS
E0    = 0.45e12   # Beam energy          [eV]
epsx_n = 2.5e-6   # Normalized emittance [m*rad]
epsy_n = 2.5e-6   # Normalized emittance [m*rad]

betax =   3.688837   # [m]
betay =  40.410416
alphax = -0.3109329  # [rad]
alphay =  3.4062036  
gammax =  0.297297   # [m]
gammay =  0.311856

dx     =  1.3454253 # Dispersion [m]
dy     =  0.0       #
dpx    =  0.1133949 # [rad]
dpy    =  0.1133949 #

rmsZ = 0.0
rmsE = 0.0

x0 = 0.1345193  #[mm]
y0 = 0.0
xp0 = 0.0113387 #[mrad]
yp0 = 0.0

#General physics parameters
mp = 938.272046e6 #proton mass, eV/c^2
e  = 1.60217657e-19 #C, electron charge
c  = 2.99792485e8   #m/s, speed of light

#Calculate parameters
# Particles should come in packets of 2
assert npart % 2 == 0 
gamma_rel = E0/mp
beta_rel  = np.sqrt(1-gamma_rel**-2)
p0 = np.sqrt((E0-mp)*(E0+mp))
# Geometrical emittance
epsx_g = epsx_n / (beta_rel*gamma_rel)
epsy_g = epsy_n / (beta_rel*gamma_rel)

def normalize_inv(x_n,xp_n,beta,alpha,emittance):
    """
    Converts normalized coordinates to non-normalized coordinates
    via equation 2.68 in 'Accelerator Physics' by Lee
    """
    emittance_sqrt = np.sqrt(emittance)
    x  = np.sqrt(beta)*x_n*emittance_sqrt
    xp = (-alpha/np.sqrt(beta)*x_n + xp_n/np.sqrt(beta))*emittance_sqrt
    
    return (x,xp)

x   = np.zeros(npart)
xp  = np.zeros(npart)
y   = np.zeros(npart)
yp  = np.zeros(npart)
z   = np.zeros(npart)
dPP = np.zeros(npart)
E   = np.zeros(npart)

i = 0
for i_amp in xrange(npart_amp):
    amp = (amp_max-amp_min)*float(i_amp)/npart_amp + amp_min
    for i_ang in xrange(npart_ang):
        ang = 2*np.pi*float(i_ang)/npart_ang + ang_off
        x_n = amp*np.cos(ang)
        xp_n = 0.0
        y_n = amp*np.sin(ang)
        yp_n = 0.0
        X,XP = normalize_inv(x_n,xp_n,betax,alphax,epsx_g)
        Y,YP = normalize_inv(y_n,yp_n,betay,alphay,epsy_g)
        x[i]  = X  + x0
        xp[i] = XP + xp0
        y[i]  = Y  + y0
        xp[i] = YP + yp0

        z[i]   = 0.0
        dPP[i] = 0.0
        E[i]   = E0
        
        i+=1


#Write to file
ofile = open("fort.13", 'w')
for i in xrange(0,npart,2):
    ofile.write(str(x  [i]  *1e3 )+"\n") #mm
    ofile.write(str(xp [i]  *1e3 )+"\n") #mrad
    ofile.write(str(y  [i]  *1e3 )+"\n") #mm
    ofile.write(str(yp [i]  *1e3 )+"\n") #mrad
    ofile.write(str(z  [i]  *1e3 )+"\n") #mm
    ofile.write(str(dPP[i]       )+"\n") #-

    ofile.write(str(x  [i+1]*1e3 )+"\n") #mm
    ofile.write(str(xp [i+1]*1e3 )+"\n") #mrad
    ofile.write(str(y  [i+1]*1e3 )+"\n") #mm
    ofile.write(str(yp [i+1]*1e3 )+"\n") #mrad
    ofile.write(str(z  [i+1]*1e3 )+"\n") #mm
    ofile.write(str(dPP[i+1]     )+"\n") #-

    ofile.write(str(E0*1e-6      )+"\n") #MeV
    ofile.write(str(E  [i]  *1e-6)+"\n") #MeV
    ofile.write(str(E  [i+1]*1e-6)+"\n") #MeV
ofile.close()

import matplotlib.pyplot as plt
plt.plot(x,y,'+')
plt.show()
