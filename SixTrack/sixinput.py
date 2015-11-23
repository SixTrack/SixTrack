#!/usr/bin/ipython

"""Read and write SixTrack input files"""

import os
import gzip
from collections import OrderedDict,namedtuple
from math import factorial

import numpy as np

clight=299792458
pi=np.pi

def getlines(fn):
    if fn.endswith('.gz'):
      fh=gzip.open(fn)
    else:
      fh = file(fn)
    for ll in fh:
      if not ll.startswith('/'):
         yield ll.strip()

def myfloat(l):
  return float(l.replace('D', 'E').replace('d', 'E'))


def bn_mad(bn_mad,n,sign):
    return sign*bn_mad*factorial(n-1)

def bn_rel(bn16,bn3,r0,d0,sign):
    out=[]
    for nn,(a,b) in enumerate(zip(bn16,bn3)):
        n=nn+1;
        sixval=d0*a*b*r0**(1-n)*10**(3*n-6)
        out.append(bn_mad( sixval,n,sign))
    return out


def readf16(fn):
    if fn.endswith('.gz'):
      fh=gzip.open(fn)
    else:
      fh = file(fn)
    out=[]
    state='label'
    for lll in fh:
        if state=='label':
            bn=[]; an=[]
            name=lll.strip()
            state='data'
        elif state=='data':
            ddd=map(myfloat,lll.split())
            if len(bn)<20:
                bn.extend(ddd)
            else:
                an.extend(ddd)
            if len(an)==20:
                state='label'
                out.append((name,bn,an))
    return out

def readf8(fn):
  if fn.endswith('.gz'):
    fh=gzip.open(fn)
  else:
    fh = file(fn)
  out=[]
  for lll in fh:
    name,rest=lll.split(None,1)
    data=map(myfloat,rest.split())
    out.append((name,data))
  return out


class Variable(object):
    def __init__(self,default,block,doc):
        self.default=default
        self.vtype=type(default)
        self.block=block
        self.doc=doc
    def __repr__(self):
        tmp="Variable(%s,%s,%s)"
        return tmp%tuple(map(repr,(self.default,self.block,self.doc)))
    def validate(self,value):
        if self.vtype is float:
           return myfloat(value)
        else:
           return self.vtype(value)




class SixTrackInput(object):
  classes=dict(
    drift=namedtuple('drift',['l']),
    mult =namedtuple('mult','knl ksl hxl hyl l rel'.split()),
    cav  =namedtuple('cav','vn f lag scav'.split()),
    align=namedtuple('align','x y tilt'),
    block=namedtuple('block','elems'),
  )
  variables=OrderedDict(
  [('title', Variable('','START','Study title')),
   ('geom',  Variable('GEOM','START',
         'FREE: lattice in fort.3; GEOM: lattice in fort.2')),
   ('printout', Variable(True,'PRIN','Print input data in the outputfile')),
   ('partnum', Variable(1000.0,'BEAM','Numer of particles in bunch')),
   ('emitnx', Variable(1000.0,'BEAM','Horizontal normalized emittance')),
   ('emitny', Variable(1000.0,'BEAM','Vertical normalized emittance')),
   ('sigz', Variable(1000.0,'BEAM','r.m.s bunch length')),
   ('sige', Variable(1000.0,'BEAM','r.m.s energy spread')),
   ('ibeco', Variable(1,'BEAM','Switch to subtract the closed orbit')),
   ('ibtyp', Variable(
          1,'BEAM','Switch to use the fast beam-beam algorithms')),
   ('lhc', Variable(0,'BEAM','Switch for anti-symmetric IR')),
   ('ibbc', Variable(0,'BEAM','Switch for linear coupling in 4D and 6D')),
   ('ctype', Variable(0,'CORR','Correction type')),
   ('ncor',  Variable(0,'CORR',
         'Number of zero length elements to be used as correctors')),
   ('comment', Variable('','COMM','Comment line')),
   ('deco_name1', Variable('','DECO','Name of skew-quadrupole family 1')),
   ('deco_name2', Variable('','DECO','Name of skew-quadrupole family 2')),
   ('deco_name3', Variable('','DECO','Name of skew-quadrupole family 3')),
   ('deco_name4', Variable('','DECO','Name of skew-quadrupole family 4')),
   ('deco_name5', Variable('','DECO','Name of focusing quadrupole families')),
   ('deco_name6', Variable('','DECO','Name of defocusing quadrupole families')),
   ('deco_Qx', Variable(0,'DECO','Horizontal tune including the integer part')),
   ('deco_Qy', Variable(0,'DECO','Vertical tune including the integer part')),
   ('diff_nord', Variable(0,'DIFF','Order of the map')),
   ('diff_nvar', Variable(0,'DIFF','Number of the variables')),
   ('diff_preda', Variable(1e-38,'DIFF','Precision needed by the DA package')),
   ('diff_nsix',
    Variable(0,'DIFF','Switch to calculate 5x6 instead of 6x6 map')),
   ('diff_ncor', Variable(0,'DIFF','Number of zero-length elements')),
   ('diff_name', Variable('','DIFF','Names of zero-length elements')),
   ('izu0', Variable(100000,'FLUC',
                     'Start value for the random number generator')),
   ('mmac', Variable(1,'FLUC','Disabled parameter, fixed to be 1')),
   ('mout', Variable(7,'FLUC','Binary switch for various purposes')),
   ('mcut', Variable(3,'FLUC','Random distribution cut of sigma')),
   ('itra', Variable(0,'INIT','Number of particles')),
   ('chi0', Variable(0.0,'INIT','Starting phase of initial coordinates')),
   ('chid',
    Variable(0.0,'INIT','Phase difference between first and second particle')),
   ('rat',
    Variable(999,'INIT','Emittance ratio of horiz. & vertical motion')),
   ('iver', Variable(0,'INIT','Switch to set vertical coordinates to zero')),
   ('itco', Variable(50,'ITER','Number of closed orbit search iterations')),
   ('dma', Variable(1e-12,'ITER','Precision of closed orbit displacement')),
   ('dmap', Variable(1e-15,'ITER','Precision of closed orbit divergence')),
   ('itqv', Variable(10,'ITER','Q adjustment iterations')),
   ('dkq', Variable(1e-10,'ITER','Q adjustment quad strength steps')),
   ('dqq', Variable(1e-10,'ITER','Q adjustment precision on tunes')),
   ('itcro', Variable(10,'ITER','Chromaticity: correction iterations')),
   ('dsm0', Variable(1e-10,'ITER','Chromaticity: sextupole strengths step')),
   ('dech', Variable(1e-10,'ITER','Chromaticity: correction precision')),
   ('de0',
    Variable(1e-09,'ITER','Variation of momentum spread for chrom')),
   ('ded',
    Variable(1e-09,'ITER','Variation of momentum spread for disp')),
   ('dsi', Variable(1e-09,'ITER','Desired orbit r.m.s value')),
   ('aper1', Variable(1000,'ITER','Horizontal aperture limit [mm]')),
   ('aper2', Variable(1000,'ITER','Vertical aperture limit [mm]')),
   ('mode',
    Variable('ELEMENT','LINE','Printout after each single ELEMENT or BLOCK')),
   ('number_of_blocks',
    Variable(0,'LINE','number of the blocks in the structure')),
   ('ilin', Variable(0,'LINE','Logical switch to calculate linear optics')),
   ('ntco',
    Variable(0,'LINE','Swtich to write out linear coupling parameters')),
   ('E_I', Variable(0,'LINE','Eigen emittance 1')),
   ('E_II', Variable(0,'LINE','Eigen emittance 2')),
   ('nord', Variable(0,'NORM','Order of the normal form')),
   ('nvar', Variable(0,'NORM','Number of varibles')),
   ('sigmax',
    Variable(0,'ORBI','Desired r.m.s for randomly distributed closed orbit')),
   ('sigmay',
    Variable(0,'ORBI','Desired r.m.s for randomly distrubuted closed orbit')),
   ('ncorru', Variable(0,'NORM','Number of correctors to be used')),
   ('ncorrep', Variable(0,'NORM','Number of corrections')),
   ('post_comment', Variable('','POST','Postprocessing comment title')),
   ('post_iav', Variable(20,'POST',
        'Averaging interval of the values of the distance in phase space')),
   ('post_nstart', Variable(0,'POST','Start turn number for the analysis')),
   ('post_nstop', Variable(0,'POST','Stop turn number for the analysis')),
   ('post_iwg', Variable(1,'POST',
       'Switch for the weighting of the slope of the distance in phase ')),
   ('post_dphix', Variable(0.08,'POST','')),
   ('post_dphiy', Variable(0.08,'POST','')),
   ('post_iskip', Variable(1,'POST','')),
   ('post_iconv', Variable(0,'POST','')),
   ('post_imad', Variable(0,'POST','')),
   ('post_cma1', Variable(1.0,'POST','')),
   ('post_cma2', Variable(1.0,'POST','')),
   ('post_Qx0', Variable(62.0,'POST','')),
   ('post_Qy0', Variable(60.0,'POST','')),
   ('post_ivox', Variable(1,'POST','')),
   ('post_ivoy', Variable(1,'POST','')),
   ('post_ires', Variable(10,'POST','')),
   ('post_dres', Variable(0.005,'POST','')),
   ('post_ifh', Variable(1,'POST','')),
   ('post_dfft', Variable(0.05,'POST','')),
   ('post_kwtype', Variable(0,'POST','')),
   ('post_itf', Variable(1,'POST','')),
   ('post_icr', Variable(0,'POST','')),
   ('post_idis', Variable(1,'POST','')),
   ('post_icow', Variable(1,'POST','')),
   ('post_istw', Variable(1,'POST','')),
   ('post_iffw', Variable(1,'POST','')),
   ('post_nprint', Variable(1,'POST','')),
   ('post_ndafi', Variable(30,'POST','')),
   ('resonance_nr', Variable(0,'RESO','')),
   ('resonance_n', Variable(0,'RESO','')),
   ('resonance_ny1', Variable(0,'RESO','')),
   ('resonance_ny2', Variable(0,'RESO','')),
   ('resonance_ny3', Variable(0,'RESO','')),
   ('resonance_ip1', Variable(0,'RESO','')),
   ('resonance_ip2', Variable(0,'RESO','')),
   ('resonance_ip3', Variable(0,'RESO','')),
   ('resonance_nrs', Variable(0,'RESO','')),
   ('resonance_ns1', Variable(0,'RESO','')),
   ('resonance_ns2', Variable(0,'RESO','')),
   ('resonance_ns3', Variable(0,'RESO','')),
   ('resonance_length', Variable(0,'RESO','')),
   ('resonance_Qx', Variable(0,'RESO','')),
   ('resonance_Qy', Variable(0,'RESO','')),
   ('resonance_Ax', Variable(0,'RESO','')),
   ('resonance_Ay', Variable(0,'RESO','')),
   ('resonance_name1', Variable('','RESO','')),
   ('resonance_name2', Variable('','RESO','')),
   ('resonance_name3', Variable('','RESO','')),
   ('resonance_name4', Variable('','RESO','')),
   ('resonance_name5', Variable('','RESO','')),
   ('resonance_name6', Variable('','RESO','')),
   ('resonance_name7', Variable('','RESO','')),
   ('resonance_name8', Variable('','RESO','')),
   ('resonance_name9', Variable('','RESO','')),
   ('resonance_name10', Variable('','RESO','')),
   ('resonance_nch', Variable(0,'RESO','')),
   ('resonance_nq', Variable(0,'RESO','')),
   ('resonance_Qx0', Variable(0,'RESO','')),
   ('resonance_Qy0', Variable(0,'RESO','')),
   ('numl',
    Variable(0,'TRAC','Number of turns in the forward direction')),
   ('numlr',
    Variable(0,'TRAC','Number of turns in the backward direction')),
   ('napx', Variable(999,'TRAC','Number of amplitude variations')),
   ('amp1',
    Variable(999,'TRAC','Start amplitude in the horiz. phase space')),
   ('amp0', Variable(999,'TRAC','End amplitude in the horiz. phase space')),
   ('ird', Variable(0,'TRAC','Switch for the type of amplitude variation')),
   ('imc',
    Variable(1,'TRAC','Number of variations of delta')),
   ('idy1', Variable(1,'TRAC','Switch for turning horiz. coupling on/off')),
   ('idy2', Variable(1,'TRAC','Switch for turning vertical coupling on/off')),
   ('idfor',
    Variable(0,'TRAC','Switch to add closed orbit to initial coordinates')),
   ('irew', Variable(1,'TRAC','Switch to save all or some tracking data')),
   ('iclo6',
    Variable(2,'TRAC','Switch to calculate 6D closed orbit using DA-package')),
   ('nde1', Variable(0,'TRAC','Number of turns at flat bottom')),
   ('nde2', Variable(0,'TRAC','Number of turns for the energy ramping')),
   ('nwr1',
    Variable(1,'TRAC',"Coordinate writing interval during flat bottom")),
   ('nwr2',
    Variable(1,'TRAC',"Coordinate writing interval during ramp")),
   ('nwr3',
    Variable(1,'TRAC',"Coordinate writing interval during flat top")),
   ('nwr4', Variable(1,'TRAC','Coordinate writing interval in fort.6')),
   ('ntwin', Variable(2,'TRAC','For analysis of Lyapunov exponent')),
   ('ibidu', Variable(1,'TRAC','Switch to create or read binary dump of accelerator')),
   ('iexact', Variable(0,'TRAC','Switch to use exact drift tracking')),
  ])
  def var_from_line(self,line,vvv):
    for val,name in zip(line.split(),vvv.split()):
      setattr(self,name,self.variables[name].validate(val))
  def __init__(self,basedir='.'):
    self.basedir=basedir
    self.filenames={}
    for n in [2,3,8,16]:
      fname='fort.%d'%n
      ffname=os.path.join(basedir,fname)
      if not os.path.isfile(ffname):
         ffname+='.gz'
         if not os.path.isfile(ffname):
           ffname=None
      self.filenames[fname]=ffname
    f3 = getlines(self.filenames['fort.3'])
    while 1:
      try:
          ll = f3.next().strip()
      except StopIteration:
          break

      # START AND END BLOCKS
      if ll.startswith('FREE') or ll.startswith('GEOM'):
          self.title = ll.split(' ',1)[1]
          self.geom = ll[:4]

      elif ll.startswith('ENDE'):
          if self.geom == 'GEOM':
              f3 = getlines(self.filenames['fort.2'])
          else:
              break

      # BLOCKS FOR FORT.3 IN ALPHABETICAL ORDER
      elif ll.startswith('BEAM'):
          ll = f3.next().strip()
          lls = ll.split()
          vvv='partnum emitnx emitny sigz sige ibeco ibtyp lhc ibbc'
          self.var_from_line(ll,vvv)
          # loop over all beam-beam elements
          ll = f3.next().strip()
          self.bbelements = {}
          while not ll.startswith('NEXT'):
              name, data = ll.split(' ', 1)
              data = data.split()
              data = [int(data[0]), float(data[1]), float(data[2])]
              self.bbelements[name.strip()] = data
              ll = f3.next().strip()

      elif ll.startswith('CHRO'):
          ll = f3.next()
          self.chromcorr = {}
          while not ll.startswith('NEXT'):
              name, data = ll.split(' ', 1)
              data = data.split()
              if len(data) == 2:
                  data = [myfloat(data[0]), int(data[1])]
              else:
                  data = [myfloat(data[0])]
              self.chromcorr[name.strip()] = data
              ll = f3.next()

      elif ll.startswith('CORR'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.ctype = int(lls[0])
              self.ncor = int(lls[1])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              self.corr_names = ll.split()
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.corr_parameters = [int(lls[0]), int(lls[1]), myfloat(lls[2]), \
                          myfloat(lls[3]), myfloat(lls[4])]

      elif ll.startswith('COMB'):
          ll = f3.next()
          self.e0 = []
          self.eRpairs = []
          while not ll.startswith('NEXT'):
              e0, eRpairs = ll.split(' ', 1)
              eRpairs = eRpairs.split()
              for ii in range(0, len(eRpairs)/2):
                  eRpairs[ii*2] = myfloat(eRpairs[ii*2])
              self.e0.append(e0)
              self.eRpairs.append(eRpairs)
              ll = f3.next()
          print self.e0
          print self.eRpairs

      elif ll.startswith('COMM'):
          self.comment = ll[:]

      elif ll.startswith('DECO'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.deco_name1 = lls[0]
              self.deco_name2 = lls[1]
              self.deco_name3 = lls[2]
              self.deco_name4 = lls[3]
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.deco_name5 = lls[0]
              self.deco_Qx = myfloat(lls[1])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.deco_name6 = lls[0]
              self.deco_Qy = myfloat(lls[1])

      elif ll.startswith('DIFF'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.diff_nord = int(lls[0])
              self.diff_nvar = int(lls[1])
              self.diff_preda = myfloat(lls[2])
              self.diff_nsix = int(lls[3])
              self.diff_ncor = int(lls[4])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              if self.diff_ncor != len(ll.split()):
                  # temporary error handling, does not abort execution
                  print 'ERROR: Mismatch between #items in row 2 and ncor'
              self.diff_name = ll.split()

      elif ll.startswith('DISP'):
          ll = f3.next()
          self.displacements = {}
          while not ll.startswith('NEXT'):
              name, data = ll.split(' ', 1)
              data = [myfloat(item) for item in data.split()]
              self.displacements[name.strip()] = data

      elif ll.startswith('FLUC'):
          ll = f3.next().strip()
          lls = ll.split()
          self.izu0 = int(lls[0])
          self.mmac = int(lls[1])
          self.mout = int(lls[2])
          self.mcut = int(lls[3])

      elif ll.startswith('INIT'):
          ll = f3.next().strip()
          self.initialconditions = []
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.itra = int(lls[0])
              self.chi0 = myfloat(lls[1])
              self.chid = myfloat(lls[2])
              self.rat  = myfloat(lls[3])
              self.iver = int(lls[4])
              ll = f3.next().strip()
          while not ll.startswith('NEXT'):
              self.initialconditions.append(myfloat(ll))
              ll = f3.next().strip()

      elif ll.startswith('ITER'):
          ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.itco = int(lls[0])
              self.dma  = myfloat(lls[1])
              self.dmap = myfloat(lls[2])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.itqv = int(lls[0])
              self.dkq  = myfloat(lls[1])
              self.dqq  = myfloat(lls[2])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.itcro = int(lls[0])
              self.dsm0  = myfloat(lls[1])
              self.dech  = myfloat(lls[2])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.de0 = myfloat(lls[0])
              self.ded = myfloat(lls[1])
              self.dsi = myfloat(lls[2])
              if len(lls) == 4:
                  self.aper1 = myfloat(lls[3])
                  self.aper2 = 'not assigned'
              if len(lls) == 5:
                  self.aper1 = myfloat(lls[3])
                  self.aper2 = myfloat(lls[4])

      elif ll.startswith('LIMI'):
          ll = f3.next()
          self.aperturelimitations = {}
          while not ll.startswith('NEXT'):
              name, data = ll.split(' ', 1)
              data = data.split()
              data = [data[0], myfloat(data[1]), myfloat(data[2])]
              self.aperturelimitations[name.strip()] = data

      elif ll.startswith('LINE'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.mode = lls[0]
              self.number_of_blocks = int(lls[1])
              self.ilin = int(lls[2])
              self.ntco = int(lls[3])
              self.E_I = myfloat(lls[4])
              self.E_II = myfloat(lls[5])
              ll = f3.next()
          self.linenames = []
          while not ll.startswith('NEXT'):
              names = ll.split()
              self.linenames = self.linenames.append(names)
              ll = f3.next()

      elif ll.startswith('MULT'):
          ll = f3.next()
          try:
              self.mult
          except AttributeError:
              # self.mult not set (first occurence of MULT)
              self.mult = {}
          name, data = ll.split(' ', 1)
          data = data.split()
          data = [myfloat(data[0]), myfloat(data[1])]
          ll = f3.next()
          while not ll.startswith('NEXT'):
              lls = ll.split()
              lls = [myfloat(item) for item in lls]
              data.append(lls)
              self.mult[name.strip()] = data
              ll = f3.next()

      elif ll.startswith('NORM'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.nord = int(lls[0])
              self.nvar = int(lls[1])

      elif ll.startswith('ORBI'):
          ll = f3.next()
          lls = ll.split()
          self.sigmax  = myfloat(lls[0])
          self.sigmay  = myfloat(lls[1])
          self.ncorru  = int(lls[2])
          self.ncorrep = int(lls[3])
          ll = f3.next()
          while not ll.startswith('NEXT'):
              # not implemented, see manual at 3.5.4 Orbit Correction
              ll = f3.next()

      elif ll.startswith('ORGA'):
          ll = f3.next()
          self.organisation_ran_numb = {}
          ii = 0
          while not ll.startswith('NEXT'):
              lls = ll.split()
              name = 'orga' + str(ii)
              self.organisation_ran_numb[name] = lls
              ll = f3.next()
              ii += 1

      elif ll.startswith('POST'):
          ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              self.post_comment = ll
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.post_iav = int(lls[0])
              self.post_nstart = int(lls[1])
              self.post_nstop = int(lls[2])
              self.post_iwg = int(lls[3])
              self.post_dphix = myfloat(lls[4])
              self.post_dphiy = myfloat(lls[5])
              self.post_iskip = int(lls[6])
              self.post_iconv = int(lls[7])
              self.post_imad = int(lls[8])
              self.post_cma1 = myfloat(lls[9])
              self.post_cma2 = myfloat(lls[10])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.post_Qx0 = myfloat(lls[0])
              self.post_Qy0 = myfloat(lls[1])
              self.post_ivox = int(lls[2])
              self.post_ivoy = int(lls[3])
              self.post_ires = int(lls[4])
              self.post_dres = myfloat(lls[5])
              self.post_ifh = int(lls[6])
              self.post_dfft = myfloat(lls[7])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.post_kwtype = int(lls[0])
              self.post_itf = int(lls[1])
              self.post_icr = int(lls[2])
              self.post_idis = int(lls[3])
              self.post_icow = int(lls[4])
              self.post_istw = int(lls[5])
              self.post_iffw = int(lls[6])
              self.post_nprint = int(lls[7])
              self.post_ndafi = int(lls[8])

      elif ll.startswith('PRIN'):
          self.printout = True
          f3.next()

      elif ll.startswith('RESO'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.resonance_nr = int(lls[0])
              self.resonance_n = int(lls[1])
              self.resonance_ny1 = int(lls[2])
              self.resonance_ny2 = int(lls[3])
              self.resonance_ny3 = int(lls[4])
              self.resonance_ip1 = int(lls[5])
              self.resonance_ip2 = int(lls[6])
              self.resonance_ip3 = int(lls[7])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.resonance_nrs = int(lls[0])
              self.resonance_ns1 = int(lls[1])
              self.resonance_ns2 = int(lls[2])
              self.resonance_ns3 = int(lls[3])
              ll = f3.next()
          if not startswith('NEXT'):
              lls = ll.split()
              self.resonance_length = myfloat(lls[0])
              self.resonance_Qx = myfloat(lls[1])
              self.resonance_Qy = myfloat(lls[2])
              self.resonance_Ax = myfloat(lls[3])
              self.resonance_Ay = myfloat(lls[4])
              ll = f3.next()
          if not startswith('NEXT'):
              lls = ll.split()
              self.resonance_name1 = lls[0]
              self.resonance_name2 = lls[1]
              self.resonance_name3 = lls[2]
              self.resonance_name4 = lls[3]
              self.resonance_name5 = lls[4]
              self.resonance_name6 = lls[5]
              ll = f3.next()
          if not startswith('NEXT'):
              lls = f3.next()
              self.resonance_nch = int(lls[0])
              self.resonance_name7 = lls[1]
              self.resonance_name8 = lls[2]
              ll = f3.next()
          if not startswith('NEXT'):
              lls = f3.next()
              self.resonance_nq = int(lls[0])
              self.resonance_name9 = lls[1]
              self.resonance_name10 = lls[2]
              self.resonance_Qx0 = myfloat(lls[3])
              self.resonance_Qy0 = myfloat(lls[4])

      elif ll.startswith('RIPP'):
          ll = f3.next()
          self.ripp = {}
          while not ll.startswith('NEXT'):
              name, settings = ll.split(' ', 1)
              settings = settings.split()
              data = [myfloat(item) for item in settings[0:3]]
              data.append(int(settings[3]))
              self.ripp[name] = data
              ll = f3.next()

      elif ll.startswith('SEAR'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.search_Qx = myfloat(lls[0])
              self.search_Qy = myfloat(lls[1])
              self.search_Ax = myfloat(lls[2])
              self.search_Ay = myfloat(lls[3])
              self.search_length = myfloat(lls[4])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.search_npos = int(lls[0])
              self.search_n = int(lls[1])
              self.search_ny1 = int(lls[2])
              self.search_ny2 = int(lls[3])
              self.search_ny3 = int(lls[4])
              self.search_ip1 = int(lls[5])
              self.search_ip2 = int(lls[6])
              self.search_ip3 = int(lls[7])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              self.sear_name = ll.split()

      elif ll.startswith('SUBR'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.subres_n1 = int(lls[0])
              self.subres_n2 = int(lls[1])
              self.subres_Qx = myfloat(lls[2])
              self.subres_Qy = myfloat(lls[3])
              self.subres_Ax = myfloat(lls[4])
              self.subres_Ay = myfloat(lls[5])
              self.subres_Ip = int(lls[6])
              self.subres_length = myfloat(lls[7])

      elif ll.startswith('SYNC'):
          ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.harm  = myfloat(lls[0])
              self.alc   = myfloat(lls[1])
              self.u0    = myfloat(lls[2])
              self.phag  = myfloat(lls[3])
              self.tlen  = myfloat(lls[4])
              self.pma   = myfloat(lls[5])
              self.ition = int(lls[6])
              if len(lls) == 8:
                  self.dppoff = myfloat(lls[7])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.dpscor = myfloat(lls[0])
              self.sigcor = myfloat(lls[1])

      elif ll.startswith('TRAC'):
          ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.numl  = int(lls[0])
              self.numlr = int(lls[1])
              self.napx  = int(lls[2])
              self.amp1  = myfloat(lls[3])
              self.amp0  = myfloat(lls[4])
              self.ird   = int(lls[5])
              self.imc   = int(lls[6])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.idy1  = int(lls[0])
              self.idy2  = int(lls[1])
              self.idfor = int(lls[2])
              self.irew  = int(lls[3])
              self.iclo6 = int(lls[4])
              ll = f3.next().strip()
          if not ll.startswith('NEXT'):
              vvv='nde1 nde2 nwr1 nwr2 nwr3 nwr4 ntwin ibidu iexact'
              self.var_from_line(ll,vvv)

      elif ll.startswith('TUNE'):
          ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.tune_name1 = lls[0]
              self.tune_Qx = myfloat(lls[1])
              if len(lls) == 3:
                  self.tune_iqmod6 = int(lls[2])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.tune_name2 = lls[0]
              self.tune_Qy = myfloat(lls[1])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.tune_name3 = lls[0]
              self.tune_deltaQ = myfloat(lls[1])
              ll = f3.next()
          if not ll.startswith('NEXT'):
              lls = ll.split()
              self.tune_name4 = lls[0]
              self.tune_name5 = lls[1]

      elif ll.startswith('TROM'):
          ll = f3.next()
          self.trom = {}
          while not ll.startswith('NEXT'):
              name = ll
              ll = f3.next()
              cvar = [myfloat(item) for item in ll.split()]
              ll = f3.next()
              cvar = cvar + [myfloat(item) for item in ll.split()]
              Mvar = []
              for ii in range(0, 12):
                  ll = f3.next()
                  Mvar = Mvar + [myfloat(item) for item in ll.split()]
              data = cvar + Mvar
              self.trom[name] = data
              ll = f3.next()

      # BLOCKS FOUND IN FORT.2 (IF GEOM)
      elif ll.startswith('SING'):
          self.single = {}
          ll = f3.next().strip()
          while not ll.startswith('NEXT'):
              name, etype, data = ll.split(None, 2)
              data = [myfloat(dd) for dd in data.split()]
              self.single[name.strip()] = [int(etype)]+data
              ll = f3.next().strip()

      elif ll.startswith('BLOC'):
          lls = f3.next().strip().split()
          self.mper = lls[0]
          self.msym = [int(ii) for ii in lls[1:]]
          ll = f3.next().strip()
          self.blocks = {}
          while not ll.startswith('NEXT'):
              lls = ll.split()
              self.blocks[lls[0]] = lls[1:]
              ll = f3.next().strip()

      elif ll.startswith('STRU'):
          self.struct = []
          ll = f3.next().strip()
          while not ll.startswith('NEXT'):
              self.struct.extend(ll.split())
              ll = f3.next().strip()
    self.add_default_vars()
    #self.add_struct_count()
    for nnn,data in self.mult.items():
        rref,bend = data[:2]
        bnrms,bn,anrms,an=zip(*data[2:])
        self.mult[nnn]=(rref,bend,bn,an,bnrms,anrms)
    if 'fort.16' in self.filenames:
       self.multblock={}
       for name,bn,an in readf16(self.filenames['fort.16']):
          self.multblock.setdefault(name,[]).append((bn,an))
    if 'fort.8' in self.filenames:
       self.align={}
       for name,(dx,dy,tilt) in readf8(self.filenames['fort.8']):
           self.align.setdefault(name,[]).append((dx,dy,tilt))
    print self.prettyprint(full=False)
  def add_default_vars(self):
      for name,var in self.variables.items():
          if not hasattr(self,name):
              setattr(self,name,var.default)
  def add_struct_count(self):
    count={}
    out=[]
    for nnn in self.struct:
        ccc=count.setdefault(nnn,0)
        out.append((nnn,ccc))
        count[nnn]=ccc+1
    self.struct=out
  def __repr__(self):
    return "<SixTrackInput %s>"%self.basedir
  def prettyprint(self,full=False):
    """Pretty print input definitions which are different from default
       values unless full=True"""
    out = []
    for name,var in self.variables.items():
      val = getattr(self,name)
      if full or var.default!= val:
         stm="%s=%s"%(name,repr(val))
         line="%-20s #[%s] %s"%(stm,var.block,var.doc)
         if len(line)>80:
            out.append('#[%s] %s'%(var.block,var.doc))
            out.append(stm)
         else:
            out.append(line)
    out.append('Data:')
    for nnn in 'single blocks struct mult align'.split():
       if hasattr(self,nnn):
           out.append("%-20s %d"%(nnn,len(getattr(self,nnn))))
    return '\n'.join(out)
  def get_knl(self,name,count):
      bnv,anv=self.multblock[name][count]
      rref,bend,bn,an,bnrms,anrms=self.mult[name]
      cstr,cref=self.single[name][1:3]
      knl=bn_rel(bnv,bn,rref,bend,-1)
      ksl=bn_rel(anv,an,rref,bend,1)
      return knl,ksl
  def compare_madmult(s,sixname,sixcount,err,madname):
      knlmad,kslmad=err.errors_kvector(np.where(err//madname)[0][0],20)
      knl,ksl=s.get_knl(sixname,sixcount)
      res=0;cc=0
      for n,(a,b) in enumerate(zip(knlmad,knl)):
         if a!=0:
           print a,b,a/b
           res+=a/b; cc+=1
      for n,(a,b) in enumerate(zip(kslmad,ksl)):
         if a!=0:
           print a,b,a/b
           res+=a/b; cc+=1
      return 1-res/cc
  def iter_struct(self):
      for el in self.struct:
        if el in self.blocks:
          for ell in self.blocks[el]:
              yield ell
        else:
          yield el
  def expand_struct(self,convert=classes):
      out=[]
      count={}
      rest=[]
      drift=convert['drift']
      mult =convert['mult']
      cav  =convert['cav']
      align=convert['align']
      block=convert['block']
      for nnn in self.iter_struct():
          ccc=count.setdefault(nnn,0)
          etype,d1,d2,d3,d4,d5,d6=self.single[nnn]
          elem=None
          if etype in [0,25]:
              elem=drift(l=d3)
          elif abs(etype) in [1,2,3,4,5,7,8,9,10]:
              bn_six=d1;nn=abs(etype); sign=-etype/nn
              madval=bn_mad(bn_six,nn,sign)
              knl=[0]*(nn-1)+[madval]; ksl=[0]*nn
              if sign==1:
                 knl,ksl=ksl,knl
              elem=mult(knl,ksl,0,0,0,0)
          elif etype==11:
              knl,ksl=self.get_knl(nnn,ccc)
              hxl=0; hyl=0;l=0
              if d3==-1:
                  hxl=-d1; l=d2
                  knl[0]=hxl
              elif d3==-2:
                  hyl=d1; l=d2
                  ksl[0]=hxl
              elem=mult(knl,ksl,hxl,hyl,l,0)
          elif etype==12:
              e0=self.initialconditions[-1]
              p0c=np.sqrt(e0**2-self.pma**2)
              beta0=p0c/e0
              v=d1; freq=d2/self.tlen*clight*beta0
              elem=cav(v/p0c,freq,lag=180-d3,scav=-1)
          else:
              rest.append([nnn]+self.single[nnn])
          if elem is not None:
            if nnn in self.align:
              dx,dy,tilt=self.align[nnn][ccc]
              tilt=tilt*180e-3/pi
              dx*=1e-3;
              dy*=1e-3;
              elem=block([align(dx,dy,-tilt),elem,align(-dx,-dy,tilt)])
            out.append([nnn,ccc,elem])
          count[nnn]=ccc+1
      return out,rest





if __name__=='__main__':
  import sys
  if len(sys.argv)==2:
     basedir=sys.argv[1]
  else:
     basedir="."
  s = SixTrackInput(basedir)
