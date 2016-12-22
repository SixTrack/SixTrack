from matplotlib.pyplot import *

r2hel1=6.928 #from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=1.1547
offsety=2.3093
thetamax=4.920e-03 # max. kick [mrad] for rr=r2hel1 from fort.3 file

figure('elens kick',figsize=(10,10))
for fnin,fnout,offx,offy in [(1,2,0,0),(2,3,offsetx,offsety),(3,4,offsetx,0),(4,5,0,offsety)]:
  subplot(2,2,fnin)
  helin=loadtxt('HEL_DUMP_%s'%fnin)
  helout=loadtxt('HEL_DUMP_%s'%fnout)
  rrin=sqrt((helin[:,3]+offx)**2+(helin[:,5]+offy)**2)
  rrout=sqrt((helout[:,3]+offx)**2+(helout[:,5]+offy)**2)
  if np.max(rrin-rrout)==0:
    plot(rrin/sig,sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2),'.',label='offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
    plot(rrin/sig,ones(len(rrin))*thetamax,'k-',label='theta_max')
    xlabel('sqrt((x+offx)**2+(y+offy)**2)/sigma ')
    ylabel('sqrt(xp**2+yp**2)=theta(r) [mrad]')
    legend(loc='best',fontsize=10)
    tight_layout()
  else:
    print 'x or y has been changed in %s - elens should only change xp,yp'%f

