import matplotlib.pyplot as plt
import numpy as np

r2hel1=6.928 # from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=-1.1547
offsety=-2.3093
theta_r2=4.920e-03 # max. kick [mrad]
names=['turn','hel1','hel2','hel3','hel4']

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
    plt.plot(rrin/sig,fff,'.',label=r'offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
    plt.plot(rrin/sig,np.ones(len(rrin))*theta_max,'k-',label=r'$\theta_{R_2}$')
    plt.plot([R2f,R2f],[0,theta_max*1.1],'g-',label=r'$n_{\mathrm{max}}$')
    plt.plot([peakT,peakT],[0,max(fff)*1.05],'r-',label=r'$n_{\mathrm{peak}}$')
    plt.xlabel(r'$n_{\sigma}=\sqrt{(x-x_{\mathrm{off}})^2+(y-y_{\mathrm{off}})^2)}$    [$\sigma$]')
    plt.ylabel(r'$\theta(r)=\sqrt{xp^2+yp^2}$ [mrad]')
    plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
    plt.grid()
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f

plt.show()
plt.close()

# plot with number of turns
dataSets={}
for name in names:
  dataSets[name]=[]
with open('dynksets.dat','r') as iFile:
  iTurn=0
  for line in iFile.readlines():
    if (line.startswith('#')): continue
    data=line.split()
    dataSets[data[1]].append(float(data[5]))
    turn=int(float(data[0])+0.001)
    if(len(dataSets['turn'])==0):
      dataSets['turn'].append(turn)
    elif (turn>dataSets['turn'][-1]):
      dataSets['turn'].append(turn)

plt.figure('elens kick',figsize=(10,10))
for name,fnin,R in zip(names[1:],[1,2,3,4],[0.5,1,1,0.5]):
  theta_max=theta_r2*R
  plt.subplot(2,2,fnin)
  plt.plot(dataSets['turn'],np.ones(len(dataSets['turn']))*theta_max,'k-',label='theta_r2')
  plt.plot(dataSets['turn'],dataSets[name],'o-',label=name)
  plt.xlabel('turn []')
  plt.ylabel(r'$\theta(r)$ [mrad]')
  plt.legend(loc='best',fontsize=10)
  plt.tight_layout()
  plt.grid()
  
  
plt.show()
plt.close()

