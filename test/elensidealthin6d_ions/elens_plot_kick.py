import matplotlib.pyplot as plt
import numpy as np

r2hel1=6.928 # from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=-1.1547
offsety=-2.3093
theta_r2=4.920e-03 # max. kick [mrad]
particlesToPlot=['p','C-12','Fe-56','Xe-129','Tl-208','Pb-207','Pb-208']
colors=['b','g','r','c','m','k','y']
# particlesToPlot=particlesToPlot[::-1]

plt.figure('elens kick',figsize=(10,10))
for fnin,fnout,offx,offy,R,R2f,peakT in [(1,2,0,0,0.5,7,7),(2,3,offsetx,offsety,1,12,10.8),(3,4,-offsetx,0,1,5,2.91604),(4,5,0,-offsety,1/2.,3,3.48995)]:
  theta_max=theta_r2*R
  plt.subplot(2,2,fnin)
  helin=np.loadtxt('HEL_DUMP_%s'%fnin)
  helout=np.loadtxt('HEL_DUMP_%s'%fnout)
  rrin=np.sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=np.sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    fff=np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2)
    for iCase in range(4):
      for iPart in range(len(particlesToPlot)):
        kk=iPart*4+iCase
        if (iCase==0):
          plt.plot(rrin[kk::28]/sig,fff[kk::28],'%s.'%(colors[iPart]),label=particlesToPlot[iPart])
        else:
          plt.plot(rrin[kk::28]/sig,fff[kk::28],'%s.'%(colors[iPart]))
    plt.plot(rrin/sig,np.ones(len(rrin))*theta_max,'k-',label=r'$\theta_{R_2}$')
    plt.plot([R2f,R2f],[0,theta_max*1.1],'g-',label=r'$n_{\mathrm{max}}$')
    plt.plot([peakT,peakT],[0,max(fff)*1.05],'r-',label=r'$n_{\mathrm{peak}}$')
    plt.xlabel(r'$n_{\sigma}=\sqrt{(x-x_{\mathrm{off}})^2+(y-y_{\mathrm{off}})^2)}$    [$\sigma$]')
    plt.ylabel(r'$\theta(r)=\sqrt{xp^2+yp^2}$ [mrad]')
    plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
    plt.grid()
    plt.title(r'offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f

plt.show()
