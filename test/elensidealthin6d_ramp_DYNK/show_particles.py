import matplotlib.pyplot as plt
import numpy as np

clight=299792458 # [m/s]

# fields in dump files
wantedFields=['x[mm]','xp[mrad]','y[mm]','yp[mrad]','sigma[mm]','(E-E0)/E0[1]']
requiredFields=wantedFields+['s[m]','particleID','turn']

# particles in the test
particlesToPlot=['p']
Zs=[1]
Ms=[9.38272310000000000E-01] # [GeV/c2]

# reference of machine
P00=450
Laccel=26658.883200 # [m]
h=35640 # []
alfa=0.000348470016216039

# reference particle
kRef=0
M0=Ms[kRef]
Z0=Zs[kRef]
P0=P00*Z0
E0=450000.9781710788E-03 # [GeV]

# what to do
lInfos=True
lOffset=False
lPlot=True
lPlotRF=True
lPng=False

# what to plot
# in the test, each particle species is replicated 4 times:
# - onMom, no betatron amplitude and with betatron amplitude;
# - offMom, no betatron amplitude and with betatron amplitude;
amplitudeClasses=['9.19324','11.49155']
momentumClasses=['6T']
firstTurnDumpName='ALL_DUMP_first'
firstElemDumpNames=['HEL_DUMP_01','HEL_DUMP_11','HEL_DUMP_21']
Mss=[]
Zss=[]
particleNames=[]
for M,Z,particle in zip(Ms,Zs,particlesToPlot):
  for momentumClass in momentumClasses:
    for amplitudeClass in amplitudeClasses:
      Mss.append(M)
      Zss.append(Z)
      particleNames.append(particle)

def pc2E(pc,M=M0,Z=Z0):
  return np.sqrt(pc**2+M0**2)

def dE2E(dE,E=E0):
  return (dE+1)*E

def E2pc(E,M=M0):
  return np.sqrt(E**2-M**2)
def E2pcPC(E,M=M0,Z=Z0):
  return E2pc(E,M=M)/Z

def betgam(Z,M,delta=0.0,P00=P00):
    return P00*Z*(1+delta)/M
def bet(betgam):
    return betgam/np.sqrt(1+betgam**2)
def gam(betgam):
    return np.sqrt(1+betgam**2)
  
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
      
def plotDumpSingleTurn( particles, IDs, title, lPng=False, symbols=['o','*','.'] ):
  '''
  plotting particle x/y/sigma vs s in a given turn
  '''
  # auto-palette
  NUM_COLORS=len(IDs)
  cm = plt.get_cmap('gist_rainbow')
                          
  # x vs s, y vs s, delta vs s
  f, axarr = plt.subplots(3,sharex=True)
  for ll in range(3):
    axarr[ll].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    mm=0
    for kk in IDs:
      axarr[ll].plot([tmpT*0.001 for tmpT in particles[kk]._t],
                          particles[kk].__dict__[wantedFields[ll*2]],'%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),markeredgewidth=0.0)
      mm+=1
    axarr[ll].grid()
    axarr[ll].legend(loc='best',fontsize=10)
    axarr[ll].set_xlabel('s [km]')
    axarr[ll].set_ylabel(wantedFields[ll*2])
  
  f.suptitle(title)
  if (lPng):
    plt.savefig(('%s_%s.png'%(title,firstTurnDumpName)).replace(' ',''))
  else:
    plt.draw()

def plotDumpSingleTurnMulti( particles, kPart, lPng=False, symbols=['o','*','.'] ):
  '''
  plotting particle x/y/sigma vs s in a given turn
  '''
  # auto-palette
  kks=[kPart,kRef]
  NUM_COLORS=len(kks)
  cm = plt.get_cmap('gist_rainbow')
  plt.rcParams.update({'font.size': 7})
                          
  # x vs s, y vs s, delta vs s
  f, axarr = plt.subplots(3,len(momentumClasses)*len(amplitudeClasses),sharex=True)
  for ii in range(len(amplitudeClasses)):
    for jj in range(len(momentumClasses)):
      iCol=jj+ii*len(momentumClasses)
      for iRow in range(3):
        axarr[iRow,iCol].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
        for mm in range(len(kks)):
          kk=kks[mm]*len(momentumClasses)*len(amplitudeClasses)+jj*len(momentumClasses)+ii
          axarr[iRow,iCol].plot( [tmpT*0.001 for tmpT in particles[kk]._t],
                                 particles[kk].__dict__[wantedFields[iRow*2]],
                                 '%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),
                                 markeredgewidth=0.0)
        axarr[iRow,iCol].grid()
        axarr[iRow,iCol].legend(loc='best',fontsize=7)
        axarr[iRow,iCol].set_xlabel('s [km]')
        axarr[iRow,iCol].set_xlim(0,Laccel*0.001)
        axarr[iRow,iCol].set_ylabel(wantedFields[iRow*2])
        if (iRow==0):
          axarr[iRow,iCol].set_title('%s - %s'%(amplitudeClasses[ii],momentumClasses[jj]))
  
  f.suptitle(particlesToPlot[kPart],fontsize=12)
  if (lPng):
    plt.savefig(('%s_%s.png'%(particlesToPlot[kPart],firstTurnDumpName)).replace(' ',''))
  else:
    plt.draw()

def plotDumpSingleElement( particles, IDs, title, lPng=False, symbols=['o','*','.'] ):
  '''
  plotting particle phase spaces: hor, ver, long at a given element
  '''
  # auto-palette
  NUM_COLORS=len(IDs)
  cm = plt.get_cmap('gist_rainbow')
                          
  # x-xp, y-yp
  f, axarr = plt.subplots(2,2)
  for ll in range(2):
    axarr[0,ll].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    mm=0
    for kk in IDs:
      axarr[0,ll].plot(particles[kk].__dict__[wantedFields[ll*2]],
                      particles[kk].__dict__[wantedFields[ll*2+1]],'%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),markeredgewidth=0.0)
    axarr[0,ll].grid()
    axarr[0,ll].legend(loc='best',fontsize=10)
    axarr[0,ll].set_xlabel(wantedFields[ll*2])
    axarr[0,ll].set_ylabel(wantedFields[ll*2+1])
  # delta vs sigma
  ll=2
  axarr[1,0].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
  mm=0
  for kk in IDs:
    Ys=E2pcPC(dE2E(np.array(particles[kk].__dict__['(E-E0)/E0[1]'])),M=Mss[kk],Z=Zss[kk])/P00-1
    axarr[1,0].plot(particles[kk].__dict__[wantedFields[ll*2]],
                       Ys,'%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),markeredgewidth=0.0)
    mm+=1
  axarr[1,0].grid()
  axarr[1,0].legend(loc='best',fontsize=10)
  axarr[1,0].set_xlabel(wantedFields[ll*2])
  axarr[1,0].set_ylabel(r'$\delta$/q []')
  
  f.suptitle(title)
  if (lPng):
    plt.savefig(('%s_%s.png'%(title,firstElemDumpName)).replace(' ',''))
  else:
    plt.draw()

def plotDumpSingleElementMulti( particles, kPart, lPng=False, symbols=['o','*','.'], title='' ):
  '''
  plotting particle phase spaces: hor, ver, long at a given element
  '''
  # auto-palette
  kks=[kPart]
  NUM_COLORS=len(kks)
  cm = plt.get_cmap('gist_rainbow')
  plt.rcParams.update({'font.size': 7})
                          
  # x vs s, y vs s, delta vs s
  f, axarr = plt.subplots(3,len(momentumClasses)*len(amplitudeClasses))
  for ii in range(len(amplitudeClasses)):
    for jj in range(len(momentumClasses)):
      iCol=jj+ii*len(momentumClasses)
      # x-xp, y-yp
      for iRow in range(2):
        axarr[iRow,iCol].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
        for mm in range(len(kks)):
          kk=kks[mm]*len(momentumClasses)*len(amplitudeClasses)+jj*len(momentumClasses)+ii
          axarr[iRow,iCol].plot(particles[kk].__dict__[wantedFields[iRow*2]],
                                particles[kk].__dict__[wantedFields[iRow*2+1]],
                                '%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),
                                markeredgewidth=0.0)
        axarr[iRow,iCol].grid()
        axarr[iRow,iCol].legend(loc='best',fontsize=7)
        axarr[iRow,iCol].set_xlabel(wantedFields[iRow*2])
        axarr[iRow,iCol].set_ylabel(wantedFields[iRow*2+1])
        if (iRow==0):
          axarr[iRow,iCol].set_title('%s - %s'%(amplitudeClasses[ii],momentumClasses[jj]))
          
      # delta vs sigma
      iRow+=1
      axarr[iRow,iCol].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
      for mm in range(len(kks)):
        kk=kks[mm]*len(momentumClasses)*len(amplitudeClasses)+jj*len(momentumClasses)+ii
        Ys=E2pcPC(dE2E(np.array(particles[kk].__dict__['(E-E0)/E0[1]'])),M=Mss[kk],Z=Zss[kk])/P00-1
        axarr[iRow,iCol].plot(particles[kk].__dict__[wantedFields[iRow*2]],Ys,
                              '%s-'%(symbols[mm]),label='ID %.0f - %s'%(kk+1,particleNames[kk]),
                              markeredgewidth=0.0)
      axarr[iRow,iCol].grid()
      axarr[iRow,iCol].legend(loc='best',fontsize=7)
      axarr[iRow,iCol].set_xlabel(wantedFields[iRow*2])
      axarr[iRow,iCol].set_ylabel(r'$\delta$/q []')
  
  f.suptitle(title,fontsize=12)
  if (lPng):
    plt.savefig(('%s_%s.png'%(particlesToPlot[kPart],firstElemDumpName)).replace(' ',''))
  else:
    plt.draw()

def plotAllLong( particles, lPng=False ):
  '''
  plotting all particles in longitudinal phase space: hor, ver, long at a given element
  '''
  # auto-palette
  NUM_COLORS=len(particlesToPlot)
  cm = plt.get_cmap('gist_rainbow')
  f, axarr = plt.subplots(1,2)
  plt.rcParams.update({'font.size': 7})
                          
  # delta vs sigma
  for ii in range(len(amplitudeClasses)):
    for jj in range(len(momentumClasses)):
      axarr[ii].set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
      for mm in range(len(particlesToPlot)):
        kk=mm*len(amplitudeClasses)*len(momentumClasses)+jj*len(amplitudeClasses)+ii
        Ys=E2pcPC(dE2E(np.array(particles[kk].__dict__['(E-E0)/E0[1]'])),M=Mss[kk],Z=Zss[kk])/P00-1
        if ( jj==0 ):
          axarr[ii].plot(particles[kk].__dict__['sigma[mm]'],Ys,
                         'o-',label=particleNames[kk],markeredgewidth=0.0)
        else:
          axarr[ii].plot(particles[kk].__dict__['sigma[mm]'],Ys,
                         'o-',markeredgewidth=0.0)
    axarr[ii].set_title(amplitudeClasses[ii])
    axarr[ii].grid()
    axarr[ii].legend(loc='best',fontsize=7)
    axarr[ii].set_xlabel('sigma[mm]')
    # axarr[ii].set_xlim(-300,300)
    axarr[ii].set_ylabel(r'$\delta$/q []')
  
  f.suptitle('Longitudinal plane - all particles')
  if (lPng):
    plt.savefig('allLong.png')
  else:
    plt.draw()

if (lPlot):
  for firstElemDumpName in firstElemDumpNames:
    # - first element, all turns
    fields=None
    particlesDumpFirstElement=[]
    with open(firstElemDumpName,'r') as iFile:
      print 'reading file %s ...'%(firstElemDumpName)
      for line in iFile.readlines():
        if (line.startswith('#')):
          if ('particleID' in line):
            fields=line.split()[1:]
            for tmpField in requiredFields:
              if ( tmpField not in fields ):
                print 'header does not contain mandatory field:',tmpField
                exit()
          continue
        if (fields is None):
          print 'no header found in file!'
          exit()
    
        # read line
        data=line.split()
        for datum,field in zip(data,fields):
          if (field=='particleID'):
            ID=int(float(datum))
            if (len(particlesDumpFirstElement)==ID-1):
              particlesDumpFirstElement.append(PARTICLE())
            particlesDumpFirstElement[ID-1].fromDumpSingle(data,fields)
    
    # plot single particles
    particlesToPlotTemp=particlesToPlot[:]
    nPlots=len(amplitudeClasses)*len(momentumClasses)
    for iPlot in range(len(particlesToPlotTemp.pop(kRef))):
      print "plotting %s - Z=%.0f, M=%.3f"%(particlesToPlot[iPlot],Zs[iPlot],Ms[iPlot])
      # - first element, all turns
      plotDumpSingleElementMulti(particlesDumpFirstElement,iPlot,lPng=lPng,title=firstElemDumpName)
  plt.show()
