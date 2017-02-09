import matplotlib.pyplot as plt
import numpy as np
import sys
import glob

r2hel1=6.928 #from fort.3 [mm]
sig=r2hel1/6 # 1 sigma beam size, hel1 between 4-6 sigma
offsetx=1.1547
offsety=2.3093
thetamax=4.920e-03 # max. kick [mrad] for rr=r2hel1 from fort.3 file

assert len(sys.argv)==2
simname = sys.argv[1]

plt.figure('elens kick',figsize=(10,10))
for fnin,fnout,offx,offy in [(1,2,0,0),(2,3,offsetx,offsety),(3,4,offsetx,0),(4,5,0,offsety)]:
  plt.subplot(2,2,fnin)
  helin_path="HEL_DUMP_%s.%s" % (fnin,simname)
  #helin_path=glob.glob("HEL_DUMP_%s.*"%fnin)[0]
  helin=np.loadtxt(helin_path)
  #helout_path=glob.glob("HEL_DUMP_%s.*"%fnout)[0]
  helout_path="HEL_DUMP_%s.%s" % (fnout,simname)
  helout=np.loadtxt(helout_path)
  
  rrin=np.sqrt((helin[:,3]+offx)**2+(helin[:,5]+offy)**2)
  rrout=np.sqrt((helout[:,3]+offx)**2+(helout[:,5]+offy)**2)
  
  if np.max(rrin-rrout)==0:
    plt.plot(rrin/sig,np.sqrt((helin[:,4]-helout[:,4])**2+(helin[:,6]-helout[:,6])**2),'.',label='offx=%2.3f sigma,offy=%2.3f sigma'%(offx/sig,offy/sig))
    plt.plot(rrin/sig,np.ones(len(rrin))*thetamax,'k-',label='theta_max')
    plt.xlabel('sqrt((x+offx)**2+(y+offy)**2)/sigma ')
    plt.ylabel('sqrt(xp**2+yp**2)=theta(r) [mrad]')
    plt.legend(loc='best',fontsize=10)
    plt.tight_layout()
  else:
    print ('x or y has been changed in %s,%s - elens should only change xp,yp'%fnin,fnout)

plt.show()
