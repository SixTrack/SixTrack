import matplotlib.pyplot as plt
import numpy as np

def parseFort6(iFileName='fort.6'):
  kicks=[]
  rrs=[]
  with open(iFileName,'r') as iFile:
    for line in iFile.readlines():
      if (line.startswith('ELENS> ELENS_KICK_FOX computed at')):
        data=line.strip().split()
        rrs.append(float(data[5]))
        kicks.append(float(data[8]))
        if (kicks[-1]!=0.0):
          kicks[-1]=-kicks[-1]
  return rrs,kicks

r2hel1=6.928 # from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
theta_r2=4.920e-03 # max. kick [mrad]

oFile=open('kicks.dat','w')
rrs,kicks=parseFort6()

plt.figure('elens kick',figsize=(13,13))
for fnin,fnout,offx,offy,R,R2f,peakT in [(1,2,sig,0,0.5,7,7),(2,3,-4*sig,4*sig,1,12,10.8),(3,4,5*sig,0,1,5,2.91604),(4,5,6*sig,6*sig,1/2.61,7,3.1704)]:
  theta_max=theta_r2*R
  plt.subplot(2,2,fnin)
  helin=np.loadtxt('HEL_DUMP_%s'%fnin)
  helout=np.loadtxt('HEL_DUMP_%s'%fnout)
  rrin=np.sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=np.sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    fff=np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2)
    plt.plot(rrin,fff,'.',label=r'offx=%2.1f $\sigma$,offy=%2.1f $\sigma$'%(offx/sig,offy/sig))
    plt.plot(rrin,np.ones(len(rrin))*theta_max,'k-',label=r'$\theta_{R_2}$')
    plt.plot([R2f*sig,R2f*sig],[0,theta_max*1.1],'g-',label=r'$R_2$')
    plt.plot([peakT*sig,peakT*sig],[0,max(fff)*1.05],'r-',label=r'$n_{\mathrm{peak}}$')
    if (len(rrs)>0):
      plt.plot(rrs[fnin-1::4],kicks[fnin-1::4],'rs',label='FOX')
    plt.xlabel('r [mm]')
    # plt.xlabel(r'$n_{\sigma}=\sqrt{(x-x_{\mathrm{off}})^2+(y-y_{\mathrm{off}})^2)}$    [$\sigma$]')
    plt.ylabel(r'$\theta(r)=\sqrt{xp^2+yp^2}$ [mrad]')
    plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
    plt.grid()
    oFile.write('# %i %i \n'%(fnin,fnout))
    for tmpR,tmpF in zip(rrin,fff):
      oFile.write(' % 22.15E % 22.15E % 22.15E \n'%(tmpR,tmpR/sig,tmpF))
    oFile.write('\n\n')
  else:
    print 'x or y has been changed in %s / %s - elens should only change xp,yp'%('HEL_DUMP_%s'%fnin,'HEL_DUMP_%s'%fnout)

oFile.close()
plt.show()
