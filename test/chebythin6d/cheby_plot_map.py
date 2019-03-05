'''
color code of points:
- black: outside R (check);
- colored: inside R (check);
- magenta: outside R (check) but wronlgy labelled as inside;
- cyan: inside R (check) but wronlgy labelled as outside;
'''

import matplotlib.pyplot as plt
import numpy as np

chebyNames=['cheby1','cheby2','cheby3','cheby4']
coords=['local','map']
offx=[0,-2, 2,0]
offy=[0, 2,-2,0]
R=6.4 # [mm]
nRows=nCols=round(len(chebyNames)/2.)

pc=450E3 # [MeV]
clight=2.99792458E8   # [m/s]
mass=0.938272310E3    # [MeV/c2]
Brho=pc/(clight*1E-6) # [Tm]
betaRel=np.sqrt(1-1/(1+(pc/mass)**2))

plt.figure('cheby_map',figsize=(20,10))
for jj in range(len(chebyNames)):
  mapIn = np.loadtxt('cheby%s_pot.dat'%(jj+1))
  ids_in=np.where(mapIn[:,5]!=0)[0]   # flagged as in by 6T
  ids_out=np.where(mapIn[:,5]==0)[0]  # flagged as out by 6T
  rr=np.sqrt((mapIn[:,0]-offx[jj])**2+(mapIn[:,1]-offy[jj])**2)
  idc_in=np.where(rr<R)[0]    # actually in
  idc_out=np.where(rr>=R)[0]  # actually out

  for ii in range(len(coords)):
    plt.subplot(nRows,nCols*len(coords),ii+jj*len(coords)+1)
    plt.scatter(mapIn[:,0+ii*2][idc_out],mapIn[:,1+ii*2][idc_out],c='k',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
    if ( len(idc_out)!=0 ):
      if ( np.max(mapIn[:,5][idc_out])!=0 ):
        print ' in map of lens %s there are some points wrongly identified as belonging to domain defined by ref radius %g when they are not... '%(chebyNames[jj],R)
        idrr=np.where(mapIn[:,5][idc_out]!=0)[0]
        plt.scatter(mapIn[:,0+ii*2][idc_out][idrr],mapIn[:,1+ii*2][idc_out][idrr],c='m',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
    if ( len(idc_in)!=0 ):
      if ( np.max(mapIn[:,5][idc_in])==0 ):
        print ' in map of lens %s there are some points wrongly identified as outside domain defined by ref radius %g when they are not... '%(chebyNames[jj],R)
        idrr=np.where(mapIn[:,5][idc_in]!=0)[0]
        plt.scatter(mapIn[:,0+ii*2][idc_in][idrr],mapIn[:,1+ii*2][idc_in][idrr],c='c',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
    plt.scatter(mapIn[:,0+ii*2][idc_in],mapIn[:,1+ii*2][idc_in],c=mapIn[:,4][idc_in],edgecolors='none')#, vmin=-3E-11, vmax=3E11)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.axis('equal')
    # plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
    plt.colorbar()
    plt.grid()
    plt.clim(-191,-136)
    plt.title('%s - %s ref sys [V m]'%(chebyNames[jj],coords[ii]))
  
plt.show()
