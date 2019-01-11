import matplotlib.pyplot as plt
import numpy as np

sig=4.82065 # 1 sigma beam size, hel1 between 4-6 sigma
mm2cm=0.1

profIN=np.loadtxt('CHG1b_170523_8p75A_2-4-2kG_500V_75mA_hires_j-vs-r.txt')
f, axarr = plt.subplots(2,2)
#
axarr[0,0].grid()
axarr[0,0].set_xlabel(r'$r$ [mm]')
axarr[0,0].set_ylabel(r'$J$ [A cm$^{-2}$]')
axarr[0,0].xaxis.set_ticks(np.arange(0,90,10))
axarr[0,0].plot(profIN[:,0],profIN[:,1],'b.-')
#
axarr[0,1].grid()
axarr[0,1].set_xlabel(r'$r$ [$\sigma$]')
axarr[0,1].set_ylabel(r'$J$ [A cm$^{-2}$]')
axarr[0,1].xaxis.set_ticks(np.arange(0,20,2))
axarr[0,1].plot(profIN[:,0]/sig,profIN[:,1],'b.-')

# cumulative
cumulIN=[]
tot=0.0
for ii in range(len(profIN[:,1])):
    if (ii==0):
        rMin=0.0
    else:
        rMin=profIN[ii-1,0]
    tot+=profIN[ii,1]*np.pi*(profIN[ii,0]-rMin)*(profIN[ii,0]+rMin)*mm2cm
    cumulIN.append(tot)
cumulIN=np.array(cumulIN)
axarr[1,0].grid()
axarr[1,0].set_xlabel(r'$r$ [$\sigma$]')
axarr[1,0].set_ylabel(r'$I$ [A]')
axarr[1,0].xaxis.set_ticks(np.arange(0,20,2))
axarr[1,0].plot(profIN[:,0]/sig,cumulIN,'b.-')
# natural kick
axarr[1,1].grid()
axarr[1,1].set_xlabel(r'$r$ [$\sigma$]')
axarr[1,1].set_ylabel(r'$\theta$ [A/mm]')
axarr[1,1].xaxis.set_ticks(np.arange(0,20,2))
axarr[1,1].plot(profIN[:,0]/sig,cumulIN/profIN[:,0],'b.-')

plt.show()
