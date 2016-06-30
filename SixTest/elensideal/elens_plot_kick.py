from matplotlib.pyplot import *

helin=loadtxt('HEL1_DUMP_1')
helout=loadtxt('HEL1_DUMP_2')

r2hel1=6.928 #from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
thetamax=4.920e-03 # max. kick [mrad] for rr=r2hel1 from fort.3 file

rrin=sqrt(helin[:,3]**2+helin[:,5]**2)
rrout=sqrt(helout[:,3]**2+helout[:,5]**2)
if np.max(rrin-rrout)==0:
  plot(rrin/sig,sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2),'b.',label='e-lens kick')
  plot(rrin/sig,ones(len(rrin))*thetamax,'k-',label='theta_max')
  xlabel('sqrt(x**2+y**2)/sigma ')
  ylabel('sqrt(xp**2+yp**2)=theta(r) [mrad]')
else:
  print 'x or y has been changed in %s - elens should only change xp,yp'%f

