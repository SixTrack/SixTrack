'''
color code of points:
- black: outside R2 or inside R1 (check);
- colored: inside R2 and outside R1 (check);
- magenta: outside R2 or inside R1 (check) but wronlgy assigned a non-NULL kick;
- white: inside R2 and outside R1 (check) but wronlgy assigned a NULL kick;
'''

import matplotlib.pyplot as plt
import numpy as np

chebyNames=['cheby1','cheby2','cheby3','cheby4']
offx=[0,-2, 2,0] # [mm]
offy=[0, 2,-2,0] # [mm]
R1=[1.0,0.7,1.5,0.4] # [mm]
R2=[6.4,6.4,6.0,6.4] #[mm]
cAngles=[0,0,-90,135] # [deg]
kicks=['kx','ky','kr']
nRows=len(chebyNames)/2
nCols=len(chebyNames)/2*len(kicks)
angle=135 # [-180:180] [deg]
dAngle=10 # [deg]
particlesToPlot=['p','C-12','Fe-56','Xe-129','Tl-208','Pb-207','Pb-208']
colors=['b','g','r','c','m','k','y']

pc=450E3 # [MeV]
clight=2.99792458E8   # [m/s]
mass=0.938272310E3    # [MeV/c2]
Brho=pc/(clight*1E-6) # [Tm]
betaRel=np.sqrt(1-1/(1+(pc/mass)**2))

for label in ['nrad','kV']:

  plt.figure('cheby_kick_%s'%(label),figsize=(30,10))
  for jj in range(len(chebyNames)):
    fnin=jj+1
    fnout=jj+2
    chebin=np.loadtxt('CHEBY_DUMP_%s'%fnin)
    chebout=np.loadtxt('CHEBY_DUMP_%s'%fnout)
  
    if (np.max(np.abs(chebin[:,3]-chebout[:,3]))==0) and (np.max(np.abs(chebin[:,5]-chebout[:,5]))==0):
      rr=np.sqrt((chebin[:,3]-offx[jj])**2+(chebin[:,5]-offy[jj])**2)
      id_in=np.where(np.logical_and(rr>=R1[jj],rr< R2[jj]))[0]
      id_out=np.where(np.logical_or(rr< R1[jj],rr>=R2[jj]))[0]
      for ii in range(len(kicks)):
        plt.subplot(nRows,nCols,ii+jj*len(kicks)+1)
        if ( kicks[ii]=='kr' ):
          # radial kick:
          z=np.sqrt((chebout[:,4]-chebin[:,4])**2+(chebout[:,6]-chebin[:,6])**2)
        else:
          z=np.array(chebout[:,4+ii*2]-chebin[:,4+ii*2])
        if (label=='nrad'):
          # show nrad
          z=z*1e+6
        else:
          # show kV from mrad
          z=z*1e-3*(betaRel*clight*Brho)*1e-3
        plt.scatter(chebin[:,3][id_out],chebin[:,5][id_out],c='k'   ,edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        idr=np.where(z[id_out]!=0.0)[0]
        if ( len(idr)>0 ):
          print ' some %s kicks from lens %s outside the domain of chebyshev lens [%g,%g) are non-zero... '%(kicks[ii],chebyNames[jj],R1[jj],R2[jj])
          plt.scatter(chebin[:,3][id_out][idr],chebin[:,5][id_out][idr],c='m',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        idr=np.where(z[id_in]==0.0)[0]
        if ( len(idr)>0 ):
          print ' some %s kicks from lens %s inside the domain of chebyshev lens [%g,%g) are zero... '%(kicks[ii],chebyNames[jj],R1[jj],R2[jj])
          plt.scatter(chebin[:,3][id_in][idr],chebin[:,5][id_in][idr],c='w',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        idr=np.where(z[id_in]!=0.0)[0]
        plt.scatter(chebin[:,3][id_in][idr],chebin[:,5][id_in][idr],c=z[id_in][idr],edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        plt.axis('equal')
        # plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.colorbar()
        plt.grid()
        if (not label=='nrad'):
          if ( kicks[ii]!='kr' ):
            plt.clim(-7,7)
        else:
          if ( kicks[ii]!='kr' ):
            plt.clim(-20,20)
        plt.title('%s - %s [%s]'%(chebyNames[jj],kicks[ii],label))
    else:
      print 'x or y has been changed in %s / %s - Chebyshev lens should only change xp,yp'%('CHEBY_DUMP_%s'%fnin,'CHEBY_DUMP_%s'%fnout)
    
  plt.show()
  plt.close()

  # plot at a given angle
  plt.figure('cheby_kick_%s_angle%s'%(label,angle),figsize=(30,10))
  angMin=np.deg2rad(angle-dAngle)
  angMax=np.deg2rad(angle+dAngle)
  for jj in range(len(chebyNames)):
    fnin=jj+1
    fnout=jj+2
    chebin=np.loadtxt('CHEBY_DUMP_%s'%fnin)
    chebout=np.loadtxt('CHEBY_DUMP_%s'%fnout)
  
    if (np.max(np.abs(chebin[:,3]-chebout[:,3]))==0) and (np.max(np.abs(chebin[:,5]-chebout[:,5]))==0):
      angles=np.arctan2(chebin[:,5]-offy[jj],chebin[:,3]-offx[jj])-np.deg2rad(cAngles[jj])
      # take care of range of angles
      angles[angles<-np.pi] += 2*np.pi
      angles[angles> np.pi] -= 2*np.pi
      # angles=(2*np.pi + angles) * (angles < 0) + angles*(angles > 0)
      ids=np.where(np.logical_and(angMin<=angles,angles<=angMax))[0]
      rr=np.sqrt((chebin[:,3][ids]-offx[jj])**2+(chebin[:,5][ids]-offy[jj])**2)
      for ii in range(len(kicks)):
        plt.subplot(nRows,nCols,ii+jj*len(kicks)+1)
        if ( kicks[ii]=='kr' ):
          # radial kick:
          z=np.sqrt((chebout[:,4][ids]-chebin[:,4][ids])**2+(chebout[:,6][ids]-chebin[:,6][ids])**2)
        elif (kicks[ii]=='kx'):
          z= np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4][ids]-chebin[:,4][ids])+np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6][ids]-chebin[:,6][ids])
        elif (kicks[ii]=='ky'):
          z=-np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4][ids]-chebin[:,4][ids])+np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6][ids]-chebin[:,6][ids])
        if (label=='nrad'):
          # show nrad
          z=z*1e+6
        else:
          # show kV from mrad
          z=z*1e-3*(betaRel*clight*Brho)*1e-3
        for iCase in range(4):
          for iPart in range(len(particlesToPlot)):
            kk=iPart*4+iCase
            if (iCase==0):
              plt.plot(rr[kk::28],z[kk::28],'%so'%(colors[iPart]),label=particlesToPlot[iPart],markeredgecolor='none')
            else:
              plt.plot(rr[kk::28],z[kk::28],'%so'%(colors[iPart]),markeredgecolor='none')
        # plt.plot(rr,z,'bo')
        plt.xlabel('r [mm]')
        plt.ylabel('%s map [%s]'%(kicks[ii],label))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()
        plt.title('%s - angle=%g'%(chebyNames[jj],angle))
    else:
      print 'x or y has been changed in %s / %s - Chebyshev lens should only change xp,yp'%('CHEBY_DUMP_%s'%fnin,'CHEBY_DUMP_%s'%fnout)
    
  plt.show()
  plt.close()
