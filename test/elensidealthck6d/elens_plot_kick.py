from matplotlib.pyplot import *
from numpy import * 

# calculation of initial tracking amplitudes
betx1=19.818526 # from LINE in SixTrack
betx2=0
rat = 0.333332 # angle 15, note sqrt(eps2/eps1)=rat
phi = arctan(sqrt(rat))
eps0 = 1.5634291429807567 # [mum] 3.75 mum normalized emittance, 1.5 GeV 
eps1=eps0*(cos(phi))**2
eps2=rat**2*eps1
sigx = sqrt(betx1*eps1)+sqrt(betx2*rat*eps2)
print 'distribute particles from amp0=0 to amp1=10*sigx=%4.2f mm'%(10*sigx)

r2hel1= 28.9239 # 6*sigx, hardcoded in fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=-1.1547
offsety=-2.3093
theta_r2=4.920e-03 # max. kick [mrad] for rr=r2hel1 from fort.3 file

close('all')
figure('elens kick',figsize=(10,10))
for fnin,fnout,offx,offy in [(1,2,0,0),(2,3,offsetx,offsety),(3,4,offsetx,0),(4,5,0,offsety)]:
  subplot(2,2,fnin)
  helin=loadtxt('HEL_DUMP_%s'%fnin)
  helout=loadtxt('HEL_DUMP_%s'%fnout)
  rrin=sqrt((helin[:,3]-offx)**2+(helin[:,5]-offy)**2)
  rrout=sqrt((helout[:,3]-offx)**2+(helout[:,5]-offy)**2)
  if np.max(rrin-rrout)==0:
    plot(rrin/sig,sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2),'.',label='offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
    plot(rrin/sig,ones(len(rrin))*theta_r2,'k-',label='theta_r2')
    xlabel('sqrt((x-offx)**2+(y-offy)**2)/sigma ')
    ylabel('sqrt(xp**2+yp**2)=theta(r) [mrad]')
    legend(loc='best',fontsize=10)
    tight_layout()
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f

draw()
show()
