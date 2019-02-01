import matplotlib.pyplot as plt
import numpy as np

r2hel1=6.928 # from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=-1.1547
offsety=-2.3093
theta_r2=4.920e-04 # max. kick [mrad]
names=['turn','hel3','hel4','GLOBAL-VARS']
milli2micro=1000

plt.figure('elens kick',figsize=(10,10))
for fnin,fnout,offx,offy,R,R2f,peakT in [(1,2,0,0,0.5,7,7),(2,3,offsetx,offsety,1,12,10.8),(3,4,-offsetx,0,1,5,2.91604),(4,5,0,-offsety,1/2.,3,3.48995)]:
  theta_max=theta_r2*R
  plt.subplot(2,2,fnin)
  # first turns
  helin=np.loadtxt('HEL_DUMP_0%s'%fnin)
  helout=np.loadtxt('HEL_DUMP_0%s'%fnout)
  rrin=np.sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=np.sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    fff=np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2)
    plt.plot(rrin/sig,fff*milli2micro,'b.',label='1-2000')
    plt.xlabel(r'$n_{\sigma}=\sqrt{(x-x_{\mathrm{off}})^2+(y-y_{\mathrm{off}})^2)}$    [$\sigma_{\mathrm{inj}}$]')
    plt.ylabel(r'$\theta(r)=\sqrt{xp^2+yp^2}$ [$\mu$rad]')
    plt.plot(rrin/sig,np.ones(len(rrin))*theta_max*milli2micro,'k-')
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f
    
  # middle turns
  helin=np.loadtxt('HEL_DUMP_1%s'%fnin)
  helout=np.loadtxt('HEL_DUMP_1%s'%fnout)
  rrin=np.sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=np.sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    fff=np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2)
    plt.plot(rrin/sig,fff*milli2micro,'c.',label='50001-52000')
    plt.plot(rrin/sig,np.ones(len(rrin))*theta_max*milli2micro*0.95,'k-')
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f
    
  # last turns
  helin=np.loadtxt('HEL_DUMP_2%s'%fnin)
  helout=np.loadtxt('HEL_DUMP_2%s'%fnout)
  rrin=np.sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=np.sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    fff=np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2)
    plt.plot(rrin/sig,fff*milli2micro,'y.',label='95001-97000')
    plt.plot([R2f,R2f],[0,theta_max*milli2micro*1.1],'g-',label=r'$R_2$')
    plt.plot([peakT,peakT],[0,max(fff*milli2micro)*1.05],'r-',label=r'$n_{\mathrm{peak}}$')
    plt.plot(rrin/sig,np.ones(len(rrin))*theta_max*milli2micro*0.9,'k-',label=r'$\theta_{R_2}$')
    plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
    plt.title(r'offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
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
for name,fnin,peakVal,label in zip(names[1:],[1,2,3],[theta_r2,0.462226,450000.9781710788],[r'$\theta(r)$ [mrad]','I [A]','E [MeV]']):
  plt.subplot(2,2,fnin)
  plt.plot(dataSets['turn'],np.ones(len(dataSets['turn']))*peakVal,'k-',label='theta_r2')
  plt.plot(dataSets['turn'],dataSets[name],'o-',label=name,markeredgewidth=0.0)
  plt.xlabel('turn []')
  plt.ylabel(label)
  plt.legend(loc='best',fontsize=10)
  plt.tight_layout()
  plt.grid()
  
  
plt.show()
plt.close()

