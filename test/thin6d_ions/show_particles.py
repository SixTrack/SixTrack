import matplotlib.pyplot as plt
import numpy as np
import sys

wantedFields=['x[mm]','xp[mrad]','y[mm]','yp[mrad]','sigma[mm]','(E-E0)/E0[1]']
requiredFields=wantedFields+['s[m]','particleID','turn']

class PARTICLE():
  def __init__(self):
    for field in wantedFields:
      self.__dict__[field]=[]
    self._t =[]
  def fromDumpALL( self, data, fields ):
    atLeastOne=False
    for datum,field in zip(data,fields):
      if (field in wantedFields ):
        self.__dict__[field].append(float(datum))
        atLeastOne=True
    if atLeastOne:
      self._t.append(float(data[fields.index('s[m]')]))
  def fromDumpSingle( self, data, fields ):
    atLeastOne=False
    for datum,field in zip(data,fields):
      if (field in wantedFields ):
        self.__dict__[field].append(float(datum))
        atLeastOne=True
    if atLeastOne:
      self._t.append(int(float(data[fields.index('turn')])))
      
def plotDumpALL( particles, IDs, title ):
  '''
  plotting particle x/xp/y/yp vs s in a given turn
  '''
  f, axarr = plt.subplots(2,2)
  ii=jj=0
  for ll in range(4):
    # auto-palette
    NUM_COLORS=len(particles)
    cm = plt.get_cmap('gist_rainbow')
    for kk in IDs:
      axarr[ii,jj].set_color_cycle([cm(1.*kk/NUM_COLORS) for i in range(NUM_COLORS)])
      axarr[ii,jj].plot([tmpT*0.001 for tmpT in particles[kk]._t],
                        particles[kk].__dict__[wantedFields[ii]],'o-',label='ID %.0f'%(kk+1))
    axarr[ii,jj].grid()
    axarr[ii,jj].legend(loc='best',fontsize=10)
    axarr[ii,jj].set_xlabel('s [km]')
    axarr[ii,jj].set_ylabel(wantedFields[ll])
    ii+=1
    if ( ii%2==0 ):
      jj+=1
      ii=0
  f.suptitle(title)
  plt.draw()

def plotDumpSingle( particles, IDs, title ):
  '''
  plotting particle phase spaces at a given turn
  '''
  f, axarr = plt.subplots(3,2)
  for ii in range(3):
    # auto-palette
    NUM_COLORS=len(particles)
    cm = plt.get_cmap('gist_rainbow')
    jj=0
    for kk in IDs:
      axarr[ii,jj].set_color_cycle([cm(1.*kk/NUM_COLORS) for i in range(NUM_COLORS)])
      axarr[ii,jj].plot(particles[kk]._t,
                        particles[kk].__dict__[wantedFields[ii*2]],'o-',label='ID %.0f'%(kk+1))
    axarr[ii,jj].grid()
    axarr[ii,jj].legend(loc='best',fontsize=10)
    axarr[ii,jj].set_xlabel('turn []')
    axarr[ii,jj].set_ylabel(wantedFields[ii*2])
    jj=1
    for kk in IDs:
      axarr[ii,jj].set_color_cycle([cm(1.*kk/NUM_COLORS) for i in range(NUM_COLORS)])
      axarr[ii,jj].plot(particles[kk].__dict__[wantedFields[ii*2]],
                        particles[kk].__dict__[wantedFields[ii*2+1]],'o-',label='ID %.0f'%(kk+1))
    axarr[ii,jj].grid()
    axarr[ii,jj].legend(loc='best',fontsize=10)
    axarr[ii,jj].set_xlabel(wantedFields[ii*2])
    axarr[ii,jj].set_ylabel(wantedFields[ii*2+1])
  f.suptitle(title)
  plt.draw()

# parse files
fields=None
particles=[]
with open('ALL_DUMP_first','r') as iFile:
  print 'reading file ALL_DUMP_first ...'
  for line in iFile.readlines():
    if (line.startswith('#')):
      if ('particleID' in line):
        fields=line.split()[1:]
        for tmpField in requiredFields:
          if ( tmpField not in fields ):
            print 'header does not contain mandatory field:',tmpField
            sys.exit()
      continue
    if (fields is None):
      print 'no header found in file!'
      sys.exit()

    # read line
    data=line.split()
    for datum,field in zip(data,fields):
      if (field=='particleID'):
        ID=int(float(datum))
        if (len(particles)==ID-1):
          particles.append(PARTICLE())
        particles[ID-1].fromDumpALL(data,fields)
#
print 'plotting...'
plotDumpALL( particles, range(0,18,2), 'small amplitude' )
plotDumpALL( particles, range(1,18,2), 'big amplitude' )
plt.show()

# parse files
fields=None
particles=[]
with open('IP3_DUMP_2','r') as iFile:
  print 'reading file IP3_DUMP_2 ...'
  for line in iFile.readlines():
    if (line.startswith('#')):
      if ('particleID' in line):
        fields=line.split()[1:]
        for tmpField in requiredFields:
          if ( tmpField not in fields ):
            print 'header does not contain mandatory field:',tmpField
            sys.exit()
      continue
    if (fields is None):
      print 'no header found in file!'
      sys.exit()

    # read line
    data=line.split()
    for datum,field in zip(data,fields):
      if (field=='particleID'):
        ID=int(float(datum))
        if (len(particles)==ID-1):
          particles.append(PARTICLE())
        particles[ID-1].fromDumpSingle(data,fields)
print 'plotting...'
plotDumpSingle( particles, range(0,18,2), 'small amplitude' )
plotDumpSingle( particles, range(1,18,2), 'big amplitude' )
plt.show()
sys.exit()
