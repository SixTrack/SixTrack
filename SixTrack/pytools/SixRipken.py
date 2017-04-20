from numpy import *
from matplotlib.pyplot import *
import os
from read_fortbin import *
import gzip

# constants
mp      = 938.272046 #MeV/c^2


def gammarel(EGeV,m0=mp):
  """returns the relativistic gamma
   input: kinetic energy E [GeV], m0 [MeV]"""
  return (EGeV*1.e3+m0)/m0

def read_dump(fn,npart):
  """reads dump files for particle coordinates
  (units: x,xp,y,yp,sig,dp/p = [mm,mrad,mm,mrad,1])
  and returns a matrix p(j,i) of particle
  amplitudes with i=particle, j=turn number.
  """
  print 'read_dump: ',fn,npart
  if os.path.isfile(fn):
    nturn = len(open(fn).readlines())/npart # number of lines/number of particles
    print '... reading in file %s'%fn
    return nturn,loadtxt(fn).reshape((nturn,npart,10))
  elif os.path.isfile(fn+'.gz'):
    nturn = len(gzip.open(fn+'.gz').readlines())/npart # number of lines/number of particles
    print '... reading in file %s'%fn
    return nturn,loadtxt(fn+'.gz').reshape((nturn,npart,10))
  else:
    print 'ERROR: file %s or %s.gz does not exist!'%(fn,fn)
    return

def read_dump_norm(fn,npart):
  """reads dump files for normalized particle coordinates
  and returns the tamatrix (=matrix of eigenvectors
  = dump_tas), the inverse tamatrix (fma_tas_inv) used 
  in SixRrack for normalistaion and a matrix p 
  with the particle coordinates 
  (units: nx,nxp,ny,nyp,nsig,ndp/p = [mm,mrad,mm,mrad,1]), 
  where p(i) are the particle
  amplitudes of particle i.
  """
  # header with co,ta and tamatrix
  nheader = 76 # number of lines for header
  ta   =np.zeros(36)
  tainv=np.zeros(36)
  clo   =np.zeros(6)
  # check if file is ziped or not
  if os.path.isfile(fn):
    nturn = (len(open(fn).readlines())-nheader)/npart
  elif os.path.isfile(fn+'.gz'):
    nturn = (len(gzip.open(fn+'.gz').readlines())-nheader)/npart
  else:
    print 'ERROR: file %s or %s.gz does not exist!'%(fn,fn)
    return
  if os.path.isfile(fn):
    p     =np.loadtxt(fn).reshape((nturn,npart,10)) 
    ff = open(fn,'rb')
  elif os.path.isfile(fn+'.gz'):
    p     =np.loadtxt(fn+'.gz').reshape((nturn,npart,10)) 
    ff = gzip.open(fn+'.gz','rb')
  else:
    print 'ERROR: file %s or %s.gz does not exist!'%fn
    return
  print '... reading in file %s'%fn
  for counter in range(nheader):
    s = ff.readline().split()
    # closed orbit
    if counter == 0:
      clo = array([ float(s[i]) for i in range(2,8) ])
    # tamatrix
    if counter >=2 and counter < 38:
      ta[counter-2] = float(s[1])
    # inverse tamatrix
    if counter >=39 and counter < 75:
      tainv[counter-39] = float(s[1])
  ta=ta.reshape((6,6))
  tainv=tainv.reshape((6,6))
  return ta,tainv,clo,p

def plot_norm_phase_space(data,partid=None):
  """Plot the normalized particle trajectories. 
  If partid=None all amplitudes are plotted, 
  if partid=n with e.g. n=3 only particle n 
  is plotted."""
  # extract only data of particle *partid*
  # if partid != None, otherwise plot
  # all particles
  if partid != None:
    data = data[:,partid:partid+1,:]
  print shape(data)
  phi=arange(101)/100.0*2*pi
  # nx npx
  subplot(2,2,1)
  radius=np.max(np.max(data[:,:,3]),np.max(data[:,:,4]))
  plot(data[:,:,3],data[:,:,4],',',label='nx-npx')
  plot(radius*cos(phi),radius*sin(phi),'k-')
  xlabel(r'$x_N [1.e-3 \sqrt{\mathrm{m}}]$')
  ylabel(r'$px_N [1.e-3 \sqrt{\mathrm{m}}]$')
  # ny npy
  subplot(2,2,2)
  radius=np.max(np.max(data[:,:,5]),np.max(data[:,:,6]))
  plot(data[:,:,5],data[:,:,6],',',label='nx-npx')
  plot(radius*cos(phi),radius*sin(phi),'k-')
  xlabel(r'$y_N [1.e-3 \sqrt{\mathrm{m}}]$')
  ylabel(r'$py_N [1.e-3 \sqrt{\mathrm{m}}]$')
  # nisg ndelta
  subplot(2,2,3)
  radius=np.max(np.max(data[:,:,7]),np.max(data[:,:,8]))
  plot(data[:,:,7],data[:,:,8],',',label='nsig-ndelta')
  plot(radius*cos(phi),radius*sin(phi),'k-')
  xlabel(r'$\sigma_N [1.e-3 \sqrt{\mathrm{m}}]$')
  ylabel(r'$(\frac{\Delta p}{p_0})_N [1.e-3 \sqrt{\mathrm{m}}]$')
  # nisg (ndelta-avg(ndelta))
  subplot(2,2,4)
  radius=np.max(np.max(data[:,:,7]),np.max((data[:,:,8]-np.average(data[:,:,8]))))
  plot(data[:,:,7],data[:,:,8]-np.average(data[:,:,8]),label='nsig-(ndelta-avg(ndelta))')
  plot(radius*cos(phi),radius*sin(phi),'k-')
  xlabel(r'$\sigma_N [1.e-3 \sqrt{\mathrm{m}}]$')
  ylabel(r'$(\frac{\Delta p}{p_0})_N-\mathrm{avg}((\frac{\Delta p}{p_0})_N) [1.e-3 \sqrt{\mathrm{m}}]$')
  tight_layout()

def _plot_fft(data,partid=None):
  """plot the fft of the particle amplitudes
  where data=self.partnorm or data=self.part
  If partid=None all amplitudes are plotted, 
  if partid=n  with e.g. n=3 only particle n 
  is plotted."""
  for s in [1,2,3]:
    subplot('31%s'%s)
    for p in range(shape(data)[1]):
      if partid==None or p==partid:
        n=[data[:,p,3+(s-1)*2],data[:,p,4+(s-1)*2]]
        f=arange(len(n[0]))/float(len(n[0]))
        if(s == 1 or s==2): # counter clock wise
          plot(f,abs(fft.fft(n[0]-1j*n[1])))
          xlim(0,0.5)
        if(s == 3): # clock wise
          plot(f,abs(fft.fft(n[0]+1j*n[1])))
          xlim(0,0.1)
    xlabel('tune')
    ylabel('amplitude')
    title('Q%s'%s)
    legend(loc='best') 

class SixRipken:
  """class to check Ripken formalism implemented in SixTrack
  Args:
  -----
  basedir : directory with output files, default: .
  loc     : name of dumpfile, default: IP3_DUMP_1

  Example:
    t=SixRipken('track1/simul/62.28_60.31/0_2/e4/15/',loc='IP3_DUMP_1')
  Attributes:
  -----------
  parttot   : total number of particles (from fort.* files)
  nturn_fma : number of turns (from loc file, e.g. IP3_DUMP_1)

  fort : files from read_fortbin
  head_fort: header files from read_fortbin
  part_fort: particle amplitues from read_fortbin
  ta_fort : t.head['tamatrix'] = tamatrix from fort.90 file
             matrix is always printed out at beginning of 
             sequence, e.g. IP3, [mm,mrad,mm,mrad,mm,1]
  ta    : tamatrix from NORM_* file, e.g. NORM_IP3_DUMP_1
          [mm,mrad,mm,mrad,mm,1]
  tainv : inverse tamatrix from NORM_* file tainv=inv(ta)
  co    : closed orbit from NORM_loc file
  pnorm : normalized particle amplitudes from NORM_* file
  part  : array of particle amplitudes part[turn,particle,coordinates]
  """
  def __init__(self,basedir='.',loc='IP3_DUMP_1'):
    self.basedir = os.path.abspath(basedir)
    # from fort.* files
    self.fort    = read_allfortbin(basedir)
    self.head_fort,self.part_fort = self.fort
    self.parttot = self.head_fort['parttot'] 
    self.ta_fort = (array(self.head_fort['tamatrix'])).reshape(6,6)
    # from 'loc' file, e.g. IP3_DUMP_1 and NORM_IP3_DUMP_1
    self.nturn_fma,self.part = read_dump(os.path.join(self.basedir,loc),self.parttot)
    self.ta,self.tainv,self.co,self.partnorm=read_dump_norm(os.path.join(self.basedir,'NORM_%s'%loc),self.parttot)
  def norm(self):
    # normalize phase space
    # !!! take copy of array, otherwise self.part is modified with pco togethere !!!
    pco=np.copy(self.part)
    for p in xrange(self.parttot):
      x,xp,y,yp,sigma,delta=(pco[:,p,3:9]).T # particle coordinates [mm,mrad,mm,mrad,mm,1]
      x0,xp0,y0,yp0,sigma0,delta0=self.co # closed orbit [mm,mrad,mm,mrad,mm,1]
    #first remove closed orbit except for delta
      for z,z0 in zip([x,xp,y,yp,sigma,delta],[x0,xp0,y0,yp0,sigma0,delta0]):
        z=z-z0
    #convert to canonical variables
      px=xp*(1+delta+delta0);py=yp*(1+delta+delta0);
      z=array([x,px,y,py,sigma,delta]).T
      pco[:,p,3:9]=array([dot(self.tainv,z[i]) for i in range(len(z))])
    # convert dp/p from sqrt(m) to 1.e-3sqrt(m)
      pco[:,p,8]=pco[:,p,8]*1.e3
    return pco
  def plot_norm_dump(self,partid=None):
    """take particle amplitudes from IP*_DUMP_* files
    and normalize them with the tainv matrix. Plot the
    trajectories. If partid=None all amplitudes are 
    plotted, if partid=n with e.g. n=3 only particle n 
    is plotted."""
    # normalize
    data=self.norm()
#    print data[0,0,:]
    # plot it
    plot_norm_phase_space(data,partid=partid)
#    print data[0,0,:]
  def plot_norm(self,partid=None):
    """take particle amplitudes from NORM_IP*_DUMP_* files
    Plot the trajectories. If partid=None all amplitudes are 
    plotted, if partid=n with e.g. n=3 only particle n 
    is plotted."""
    plot_norm_phase_space(data=self.partnorm,partid=partid)
  def plot_norm_fft(self,partid=None):
    """plot fft of normalized particle coordinates
    If partid=None all amplitudes are plotted, 
    if partid=n  with e.g. n=3 only particle n 
    is plotted."""
    _plot_fft(self.partnorm,partid=partid)
  def plot_fft(self,partid=None):
    """plot fft of physical particle coordinates
    (x,x',...), delta has to rescaled to 1.e-3.
    If partid=None all amplitudes are plotted, 
    if partid=n  with e.g. n=3 only particle n 
    is plotted."""
    data=self.part
    data[:,:,8]=data[:,:,8]*1.e-3
    _plot_fft(data,partid=partid)
