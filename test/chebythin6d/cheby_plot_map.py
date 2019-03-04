import matplotlib.pyplot as plt
import numpy as np

chebyNames=['cheby1','cheby2','cheby3','cheby4']
coords=['local','map']
offx=[0,-2, 2,0]
offy=[0, 2,-2,0]
R=6.4 #[mm]
nRows=nCols=round(len(chebyNames)/2.)

pc=450E3 # [MeV]
clight=2.99792458E8   # [m/s]
mass=0.938272310E3    # [MeV/c2]
Brho=pc/(clight*1E-6) # [Tm]
betaRel=np.sqrt(1-1/(1+(pc/mass)**2))

plt.figure('cheby_map',figsize=(20,10))
for jj in range(len(chebyNames)):
  mapIn = np.loadtxt('cheby%s_pot.dat'%(jj+1))
  ids=np.where(mapIn[:,5]!=0)[0]

  for ii in range(len(coords)):
    plt.subplot(nRows,nCols*len(coords),ii+jj*len(coords)+1)
    plt.scatter(mapIn[:,0+ii*2],mapIn[:,1+ii*2],c='k',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
    plt.scatter(mapIn[:,0+ii*2][ids],mapIn[:,1+ii*2][ids],c=mapIn[:,4][ids],edgecolors='none')#, vmin=-3E-11, vmax=3E11)
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
