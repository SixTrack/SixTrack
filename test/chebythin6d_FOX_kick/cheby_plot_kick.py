'''
color code of points:
- black: outside R2 or inside R1 (check);
- colored: inside R2 and outside R1 (check);
- magenta: outside R2 or inside R1 (check) but wronlgy assigned a non-NULL kick;
- white: inside R2 and outside R1 (check) but wronlgy assigned a NULL kick;
'''

import matplotlib.pyplot as plt
import numpy as np
import sys

def parseFort6(iFileName='fort.6'):
  Xs=[]
  Ys=[]
  kickXs=[]
  kickYs=[]
  with open(iFileName,'r') as iFile:
    for line in iFile.readlines():
      if (line.startswith('CHEBY> CHEBY_KICK_FOX computing at')):
        data=line.strip().split()
        Xs.append(float(data[5]))
        Ys.append(float(data[8]))
      if (line.startswith('CHEBY> CHEBY_KICK_FOX - computed kick:')):
        data=line.strip().split()
        kickXs.append(float(data[6]))
        kickYs.append(float(data[9]))
  return np.array(Xs),np.array(Ys),np.array(kickXs),np.array(kickYs)

chebyNames=['cheby1','cheby2','cheby3','cheby4']
offx=[0,-2, 0.7,0.5]
offy=[0, 2,-0.7,0.2]
R1=[0.25,0.7,0.0,0.0] # [mm]
R2=[10.0,8.0,9.8,10.0] #[mm]
cFact=[1,1,1,2]
cAngles=[0,0,-90,160] # [deg]
kicks=['kx','ky','kr']
nRows=len(chebyNames)/2
nCols=len(chebyNames)/2*len(kicks)
dAngle=1 # [deg]
epsilon=1E-15
gteps=1+epsilon
lteps=1-epsilon
k0s=[112.422,-56.211] # [V m]
rRef=10E-03 # [m]

pc=450E3 # [MeV]
clight=2.99792458E8   # [m/s]
mass=0.938272310E3    # [MeV/c2]
Brho=pc/(clight*1E-6) # [Tm]
betaRel=np.sqrt(1-1/(1+(pc/mass)**2))

Xs,Ys,kickXs,kickYs=parseFort6()
if (len(Xs)!=len(kickXs)):
  print 'problems at parsing fort.6 file!'
  print 'found %i points, but only %i kicks!'%(len(Xs),len(kickXs))
  sys.exit(1)
if (len(Xs)==0):
  print 'found 0 FOX echoes in fort.6 file!'
  sys.exit(1)

for lKick,label in zip([True,False],['nrad','kV']):

  plt.figure('cheby_kick_%s'%(label),figsize=(30,10))
  for jj in range(len(chebyNames)):
    fnin=jj+1
    fnout=jj+2
    chebin=np.loadtxt('CHEBY_DUMP_%s'%fnin)
    chebout=np.loadtxt('CHEBY_DUMP_%s'%fnout)
  
    if (np.max(np.abs(chebin[:,3]-chebout[:,3]))==0) and (np.max(np.abs(chebin[:,5]-chebout[:,5]))==0):
      # check coordinates against R1 and R2
      xx=chebin[:,3]-offx[jj]
      yy=chebin[:,5]-offy[jj]
      angles=np.arctan2(yy,xx)-np.deg2rad(cAngles[jj])
      rr=np.sqrt(xx**2+yy**2)
      xx=rr*np.cos(angles)
      yy=rr*np.sin(angles)
      xx=np.absolute(xx)
      yy=np.absolute(yy)
      id_out=np.where(                np.logical_or( np.logical_and(xx<R1[jj]*lteps,yy<R1[jj]*lteps), np.logical_or(xx>R2[jj]*gteps,yy>R2[jj]*gteps) )) [0]
      id_in =np.where(np.logical_not( np.logical_or( np.logical_and(xx<R1[jj]*lteps,yy<R1[jj]*lteps), np.logical_or(xx>R2[jj]*gteps,yy>R2[jj]*gteps) )))[0]
      for ii in range(len(kicks)):
        plt.subplot(nRows,nCols,ii+jj*len(kicks)+1)
        if ( kicks[ii]=='kr' ):
          # radial kick:
          z=np.sqrt((chebout[:,4]-chebin[:,4])**2+(chebout[:,6]-chebin[:,6])**2)
        else:
          z=np.array(chebout[:,4+ii*2]-chebin[:,4+ii*2])
        if (label=='nrad'):
          # show nrad
          fact=1e+6
        else:
          # show kV from mrad
          fact=1e-3*(betaRel*clight*Brho)*1e-3
        z=z*fact
        # particles outside domain (black)
        plt.scatter(chebin[:,3][id_out],chebin[:,5][id_out],c='k'   ,edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        # non-zero kicks of particles outside domain (magenta)
        idr=np.where(z[id_out]!=0.0)[0]
        if ( len(idr)>0 ):
          print ' some %s kicks from lens %s outside the domain of chebyshev lens [%g,%g] are non-zero... '%(kicks[ii],chebyNames[jj],R1[jj],R2[jj])
          plt.scatter(chebin[:,3][id_out][idr],chebin[:,5][id_out][idr],c='m',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        # zero kicks of particles inside domain (white)
        idr=np.where(z[id_in]==0.0)[0]
        if ( len(idr)>0 ):
          print ' some %s kicks from lens %s inside the domain of chebyshev lens [%g,%g] are zero... '%(kicks[ii],chebyNames[jj],R1[jj],R2[jj])
          plt.scatter(chebin[:,3][id_in][idr],chebin[:,5][id_in][idr],c='w',edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        # particles inside domain (colored)
        idr=np.where(z[id_in]!=0.0)[0]
        plt.scatter(chebin[:,3][id_in][idr],chebin[:,5][id_in][idr],c=z[id_in][idr],edgecolors='none')#, vmin=-3E-11, vmax=3E11)
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        plt.axis('equal')
        # plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.colorbar()
        plt.grid()
        plt.title('%s - %s [%s]'%(chebyNames[jj],kicks[ii],label))
    else:
      print 'x or y has been changed in %s / %s - Chebyshev lens should only change xp,yp'%('CHEBY_DUMP_%s'%fnin,'CHEBY_DUMP_%s'%fnout)
    
  plt.show()
  plt.close()
  
  # # plot at a given angle
  # plt.figure('cheby_kick_%s_at_angle'%(label),figsize=(30,10))
  # for jj in range(len(chebyNames)):
  #   fnin=jj+1
  #   fnout=jj+2
  #   chebin=np.loadtxt('CHEBY_DUMP_%s'%fnin)
  #   chebout=np.loadtxt('CHEBY_DUMP_%s'%fnout)
  # 
  #   if (len(Xs)>0):
  #     angle=np.average(np.rad2deg(np.arctan2(Ys[fnin-1::4],Xs[fnin-1::4])))
  #   else:
  #     angle=160
  #   angMin=np.deg2rad(angle-dAngle)
  #   angMax=np.deg2rad(angle+dAngle)
  #   
  #   if (np.max(np.abs(chebin[:,3]-chebout[:,3]))==0) and (np.max(np.abs(chebin[:,5]-chebout[:,5]))==0):
  # 
  #     # kicks from tracking (convert to map ref sys):
  #     angles=np.arctan2(chebin[:,5]-offy[jj],chebin[:,3]-offx[jj])-np.deg2rad(cAngles[jj])
  #     # take care of range of angles
  #     angles[angles<-np.pi] += 2*np.pi
  #     angles[angles> np.pi] -= 2*np.pi
  #     # angles=(2*np.pi + angles) * (angles < 0) + angles*(angles > 0)
  #     ids=np.where(np.logical_and(angMin<=angles,angles<=angMax))[0]
  #     rr=np.sqrt((chebin[:,3][ids]-offx[jj])**2+(chebin[:,5][ids]-offy[jj])**2)
  #     for ii in range(len(kicks)):
  #       plt.subplot(nRows,nCols,ii+jj*len(kicks)+1)
  #       if ( kicks[ii]=='kr' ):
  #         # radial kick:
  #         z=np.sqrt((chebout[:,4][ids]-chebin[:,4][ids])**2+(chebout[:,6][ids]-chebin[:,6][ids])**2)
  #       elif (kicks[ii]=='kx'):
  #         z= np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4][ids]-chebin[:,4][ids])+np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6][ids]-chebin[:,6][ids])
  #       elif (kicks[ii]=='ky'):
  #         z=-np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4][ids]-chebin[:,4][ids])+np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6][ids]-chebin[:,6][ids])
  #       if (lKick):
  #         # show nrad
  #         fact=1e+6
  #       else:
  #         # show kV from mrad
  #         fact=1e-3*(betaRel*clight*Brho)*1e-3
  #       plt.plot(rr,z*fact,'bo')
  # 
  #       # kicks from FOX
  #       if (len(Xs)>0):
  #         rrs=np.sqrt(Xs**2+Ys**2) # coorindates already in map ref sys
  #         if(kicks[ii]=='kr'):
  #           yys=np.sqrt(kickXs**2+kickYs**2)
  #           plt.plot(rrs[fnin-1::4],yys[fnin-1::4]*fact,'rs',label='FOX')
  #         elif(kicks[ii]=='kx'):
  #           tmpYs= np.cos(np.deg2rad(cAngles[jj]))*kickXs[fnin-1::4]+np.sin(np.deg2rad(cAngles[jj]))*kickYs[fnin-1::4]
  #           plt.plot(rrs[fnin-1::4],tmpYs*fact,'rs',label='FOX')
  #         elif(kicks[ii]=='ky'):
  #           tmpYs=-np.sin(np.deg2rad(cAngles[jj]))*kickXs[fnin-1::4]+np.cos(np.deg2rad(cAngles[jj]))*kickYs[fnin-1::4]
  #           plt.plot(rrs[fnin-1::4],tmpYs*fact,'rs',label='FOX')
  #           
  #       plt.xlabel('r (map ref sys) [mm]')
  #       plt.ylabel('%s (map ref sys) [%s]'%(kicks[ii],label))
  #       plt.tight_layout()
  #       plt.grid()
  #       plt.title('%s - angle=%g'%(chebyNames[jj],angle))
  #   else:
  #     print 'x or y has been changed in %s / %s - Chebyshev lens should only change xp,yp'%('CHEBY_DUMP_%s'%fnin,'CHEBY_DUMP_%s'%fnout)
  #   
  # plt.show()
  # plt.close()
  
  # check linear kicks
  rrs=np.sqrt(Xs**2+Ys**2) # FOX coorindates already in map ref sys
  yys=np.sqrt(kickXs**2+kickYs**2)
  plt.figure('cheby_kick_%s_linear'%(label),figsize=(30,10))
  for jj in range(len(chebyNames)):
    fnin=jj+1
    fnout=jj+2
    chebin=np.loadtxt('CHEBY_DUMP_%s'%fnin)
    chebout=np.loadtxt('CHEBY_DUMP_%s'%fnout)

    if (np.max(np.abs(chebin[:,3]-chebout[:,3]))==0) and (np.max(np.abs(chebin[:,5]-chebout[:,5]))==0):

      # kicks from tracking (convert to map ref sys):
      angles=np.arctan2(chebin[:,5]-offy[jj],chebin[:,3]-offx[jj])-np.deg2rad(cAngles[jj])
      # take care of range of angles
      angles[angles<-np.pi] += 2*np.pi
      angles[angles> np.pi] -= 2*np.pi
      # angles=(2*np.pi + angles) * (angles < 0) + angles*(angles > 0)
      rr=np.sqrt((chebin[:,3]-offx[jj])**2+(chebin[:,5]-offy[jj])**2)
      for ii in range(len(kicks)):
        plt.subplot(nRows,nCols,ii+jj*len(kicks)+1)
        if ( kicks[ii]=='kr' ):
          # radial kick:
          tmpXs=rr
          z=np.sqrt((chebout[:,4]-chebin[:,4])**2+(chebout[:,6]-chebin[:,6])**2)
        elif (kicks[ii]=='kx'):
          tmpXs=rr*np.cos(angles)
          z= np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4]-chebin[:,4])+np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6]-chebin[:,6])
          fitXs=np.sort(tmpXs)
          fitYs=-k0s[0]/rRef**2*4*fitXs*1e-6 # fitXs are in mm, and fitYs in V
        elif (kicks[ii]=='ky'):
          tmpXs=rr*np.sin(angles)
          z=-np.sin(np.deg2rad(cAngles[jj]))*np.array(chebout[:,4]-chebin[:,4])+np.cos(np.deg2rad(cAngles[jj]))*np.array(chebout[:,6]-chebin[:,6])
          fitXs=np.sort(tmpXs)
          fitYs=-k0s[1]/rRef**2*4*fitXs*1e-6 # fitXs are in mm, and fitYs in V
        if (lKick):
          # show nrad
          fact=1e+6
          fitYs=fitYs/(1e-3*(betaRel*clight*Brho)*1e-3)*1e+6
        else:
          # show kV from mrad
          fact=1e-3*(betaRel*clight*Brho)*1e-3
        plt.plot(tmpXs,z*fact,'bo')
        if ( kicks[ii]!='kr' ):
          plt.plot(fitXs,fitYs*cFact[jj],'m-')

        # kicks from FOX
        if (len(Xs)>0):
          if(kicks[ii]=='kr'):
            plt.plot(rrs[fnin-1::4],yys[fnin-1::4]*fact,'rs',label='FOX')
          elif(kicks[ii]=='kx'):
            tmpYs= np.cos(np.deg2rad(cAngles[jj]))*kickXs[fnin-1::4]+np.sin(np.deg2rad(cAngles[jj]))*kickYs[fnin-1::4]
            plt.plot(Xs[fnin-1::4],tmpYs*fact,'rs',label='FOX')
          elif(kicks[ii]=='ky'):
            tmpYs=-np.sin(np.deg2rad(cAngles[jj]))*kickXs[fnin-1::4]+np.cos(np.deg2rad(cAngles[jj]))*kickYs[fnin-1::4]
            plt.plot(Ys[fnin-1::4],tmpYs*fact,'rs',label='FOX')
            
        if(kicks[ii]=='kr'):
          plt.xlabel('r (map ref sys) [mm]')
        elif(kicks[ii]=='kx'):
          plt.xlabel('x (map ref sys) [mm]')
        elif(kicks[ii]=='ky'):
          plt.xlabel('y (map ref sys) [mm]')
        plt.ylabel('%s (map ref sys) [%s]'%(kicks[ii],label))
        plt.tight_layout()
        plt.grid()
        plt.title('%s'%(chebyNames[jj]))
    else:
      print 'x or y has been changed in %s / %s - Chebyshev lens should only change xp,yp'%('CHEBY_DUMP_%s'%fnin,'CHEBY_DUMP_%s'%fnout)
    
  plt.show()
  plt.close()
