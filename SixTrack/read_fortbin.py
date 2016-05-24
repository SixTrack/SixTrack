# Read binary files fort.90 fort.89 ... prdouce by Sixtrack
#
# author: R. De Maria
#
#Copyright 2014 CERN. This software is distributed under the terms of the GNU
#Lesser General Public License version 2.1, copied verbatim in the file
#``COPYING''.
#
#In applying this licence, CERN does not waive the privileges and immunities
#granted to it by virtue of its status as an Intergovernmental Organization or
#submit itself to any jurisdiction.



import struct
import os
import numpy as np

def _read(fh,fmt):
  out={}
  for line in fmt.splitlines():
    lbl,spec,desc=line.split(None,2)
    data=fh.read(struct.calcsize(spec))
    obj=struct.unpack(spec,data)
    if len(obj)==1:
      obj=obj[0]
    out[lbl]=obj
  return out

fmt_head="""\
head1      1I Fortran header
title     80s General title of the run
title2    80s Additional title
date       8s Date
time       8s Time
progname   8s Program name
partfirst  1I First particle in the file
partlast   1I Last particle in the file
parttot    1I Total number of particles
spacecode  1I Code for dimensionality of phase space (1,2,4 are hor., vert. and longitudinal respectively)
turnproj   1I Projected number of turns
qx         1d Horizontal Tune
qy         1d Vertical Tune
qs         1d Longitudinal Tune
closorb    6d Closed Orbit vector
dispvec    6d Dispersion vector
tamatrix   36d Six-dimensional map for normalisation
mess1      50d 50 additional parameter
mess2      1I ...
"""
"""
seedmax    1d Maximum number of different seeds
seednum    1d Actual seed number
seedstart  1d Starting value of the seed
turnrev    1d Number of turns in the reverse direction (IBM only)
lyapcor1   1d Correction-factor for the Lyapunov (sigma=s - v0 t)
lyapcor2   1d Correction-factor for the Lyapunov (DeltaP/P0)
turnrip    1d Start turn number for ripple prolongation
"""

fmt_part="""\
partnum  1I Particle number
partdist 1d Angular distance in phase space
x        1d x (mm)
xp       1d x'(mrad)
y        1d y (mm)
yp       1d y'(mrad)
sig      1d Path-length sigma=s - v0 t (mm)
delta    1d DeltaP/P0 (1)
energy   1d Energy (Mev)
"""

def read_fortbin(fn):
  fh=open(fn,'rb')
  header=_read(fh,fmt_head)
  partfirst=header['partfirst']
  partlast=header['partlast']
  part={}
  for i in range(partfirst,partlast+1):
    part[i]=[]
  while fh.read(4)!='':  # read(fh,'headpart 1I ...')
    turnnum= struct.unpack('I',fh.read(4))
    #read(fh,fmt_part)
    #read(fh,fmt_part)
    for i in range(partfirst,partlast+1):
      pnum1= struct.unpack('I',fh.read(4))
      orb1 = struct.unpack('8d',fh.read(64))
      part[i].append(list(orb1))
    fh.read(4) # read(fh,'headpart 1I ...')
  for i in range(partfirst,partlast+1):
    part[i]=np.array(part[i])
  return header,part

def read_allfortbin(basedir='.'):
  fn=os.path.join(basedir,'fort.90')
  head,part=read_fortbin(fn)
  for lbl,obj in head.items():
    print "%-8s <%s>"%(lbl, obj)
  npart=head['parttot']
  for i in range(1,npart/2):
    fn =os.path.join(basedir,'fort.%d'%(90-i))
    nhead,ndata=read_fortbin(fn)
  #  for lbl,obj in nhead.items():
  #    print "%-8s <%s>"%(lbl, obj)
    part.update(ndata)
  for i in sorted(part.keys()):
    print "part %d turns %d"%(i,len(part[i]))
  return head,part



if __name__=='__main__':
  import sys
  if len(sys.argv)>1:
    basedir=sys.argv[1]
  else:
    basedir='.'
  print basedir
  head,part=read_allfortbin(basedir)


  #pdist,x,xp,y,yp,sigma,delta,energy=part[1].T

  #f=np.linspace(0,1,len(x))
  #tunx=np.fft.fft(x+1j*xp)
  #tuny=np.fft.fft(y+1j*yp)

  #plot(f,abs(tunx),label='qx')
  #plot(f,abs(tuny),label='qy')

#ta=array(head['tamatrix']).reshape(6,6)
#betxI=ta[0,0]**2+ta[0,1]**2
#
#J=array([[ 0.,  1.,  0.,  0.,  0.,  0.],
#         [-1.,  0.,  0.,  0.,  0.,  0.],
#         [ 0.,  0.,  0.,  1.,  0.,  0.],
#         [ 0.,  0., -1.,  0.,  0.,  0.],
#         [ 0.,  0.,  0.,  0.,  0.,  1.],
#         [ 0.,  0.,  0.,  0., -1.,  0.]])
#
#dot(ta.T,dot(J,ta))






