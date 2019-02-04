from abc import ABCMeta, abstractmethod

# Abstract class inherited by all definition classes
class Definition:
	__metaclass__ = ABCMeta
	name= ''

	@abstractmethod
	def toString():
		pass

# Only deals with pure x-poles
class MultipoleDef(Definition):
	def __init__(self, name, order, KN, KS, TILT=0, THICK=False):
		self.name = name
		self.TILT = TILT
		self.THICK = THICK

		aux = ['0' for x in range(int(order))]
		aux[-1] = str(KN)
		self.KNL = ','.join(aux)
		aux[-1] = str(KS)
		self.KSL = ','.join(aux)

	def toString(self):
		return '{}: MULTIPOLE, TILT={}, KNL={{{}}}, KSL={{{}}};'.format(
			self.name, self.TILT, self.KNL, self.KSL)

class SolenoidDef(Definition):
	def __init__(self, name, L, KS, KSI=0):
		self.name = name
		self.L = L
		self.KS = KS
		self.KSI = KSI

	def toString(self):
		return '{}: SOLENOID, L={}, KS={}, KSI={};'.format(
			self.name, self.L, self.KS, self.KSI)

# Only deals with pure RF-x-poles
class RFMultipoleDef(Definition):
	def __init__(self, name, order, VOLT, LAG, HARMON, FREQ,
		TILT, KN, KS, PN, PS):
		self.name = name
		self.VOLT = VOLT
		self.LAG = LAG
		self.HARMON = HARMON
		self.FREQ = FREQ
		self.TILT = TILT

		aux = ['0' for x in range(int(order))]
		aux[-1] = str(KN)
		self.KNL = ','.join(aux)
		aux[-1] = str(KS)
		self.KSL = ','.join(aux)
		aux[-1] = str(PN)
		self.PNL = ','.join(aux)
		aux[-1] = str(PS)
		self.PSL = ','.join(aux)

	def toString(self):
		if self.FREQ == 0:
			return '{}: RFMULTIPOLE, VOLT={}, LAG={}, HARMON={}, TILT={}, KNL={{{}}}, KSL={{{}}}, PNL={{{}}}, PSL={{{}}};'.format(
				self.name, self.VOLT, self.LAG, self.HARMON, self.TILT, self.KNL, self.KSL, self.PNL, self.PSL)
		else:
			return '{}: RFMULTIPOLE, VOLT={}, LAG={}, FREQ={}, TILT={}, KNL={{{}}}, KSL={{{}}}, PNL={{{}}}, PSL={{{}}};'.format(
				self.name, self.VOLT, self.LAG, self.FREQ, self.TILT, self.KNL, self.KSL, self.PNL, self.PSL)

class CrabCavityDef(Definition):
	def __init__(self, name, L, TILT, VOLT, LAG,
		HARMON, FREQ, RV1, RV2, RV3, RV4, LAGF, RPH1, RPH2):
		self.name = name
		self.L = L
		self.TILT = TILT
		self.VOLT = VOLT
		self.LAG = LAG
		self.HARMON = HARMON
		self.FREQ = FREQ
		self.RV1 = RV1
		self.RV2 = RV2
		self.RV3 = RV3
		self.RV4 = RV4
		self.LAGF = LAGF
		self.RPH1 = RPH1
		self.RPH2 = RPH2
		self.LAGF = LAGF

	def toString(self):
		if self.FREQ == 0:
			return '{}: CRABCAVITY, L={}, TILT={}, VOLT={}, LAG={}, HARMON={}, RV1={}, RV2={}, RV3={}, RV4={}, RPH1={}, RPH2={}, LAGF={};'.format(
				self.name, self.L, self.TILT, self.VOLT, self.LAG, self.HARMON, self.RV1, self.RV2, self.RV3, self.RV4, self.RPH1, self.RPH2, self.LAGF)
		else:
			return '{}: CRABCAVITY, L={}, TILT={}, VOLT={}, LAG={}, FREQ={}, RV1={}, RV2={}, RV3={}, RV4={}, RPH1={}, RPH2={}, LAGF={};'.format(
				self.name, self.L, self.TILT, self.VOLT, self.LAG, self.FREQ, self.RV1, self.RV2, self.RV3, self.RV4, self.RPH1, self.RPH2, self.LAGF)

class BeamBeamDef(Definition):
	def __init__(self, name, CHARGE, XMA, YMA, SIGX, SIGY, WIDTH, BBSHAPE, BBDIR):
		self.name = name
		self.CHARGE = CHARGE
		self.XMA = XMA
		self.YMA = YMA
		self.SIGX = SIGX
		self.SIGY = SIGY
		self.WIDTH = WIDTH
		self.BBSHAPE = BBSHAPE
		self.BBDIR = BBDIR

	def toString(self):
		return '{}: BEAMBEAM, CHARGE={}, XMA={}, YMA={}, SIGX={}, SIGY={}, WIDTH={}, BBSHAPE={}, BBDIR={};'.format(
			self.name, self.CHARGE, self.XMA, self.YMA, self.SIGX, self.SIGY, self.WIDTH, self.BBSHAPE, self.BBDIR)

# Could be improved upon (non-fixed L and TILT)
class KickerDef(Definition):
	def __init__(self, name, HKICK, VKICK, L, TILT):
		self.name = name
		self.HKICK = HKICK
		self.VKICK = VKICK
		self.TILT = TILT
		self.L = L

	def toString(self):
		if self.HKICK == 0:
			return '{}: VKICKER, KICK={}, L={}, TILT={};'.format(self.name,
				self.VKICK, self.L, self.TILT)
		elif self.VKICK == 0:
			return '{}: HKICKER, KICK={}, L={}, TILT={};'.format(self.name,
				self.HKICK, self.L, self.TILT)
		else:
			return '{}: KICKER, HKICK={}, VKICK={}, L={}, TILT={};'.format(
				self.name, self.HKICK, self.VKICK, self.L, self.TILT)

class RFCavityDef(Definition):
	def __init__(self, name, VOLT, LAG, L, HARMON, FREQ):
		self.name = name
		self.VOLT = VOLT
		self.LAG = LAG
		self.L = L
		self.HARMON = HARMON
		self.FREQ = FREQ

	def toString(self):
		if self.FREQ == 0:
			return '{}: RFCAVITY, VOLT={}, LAG={}, L={}, HARMON={};'.format(
				self.name, self.VOLT, self.LAG, self.L, self.HARMON)
		else:
			return '{}: RFCAVITY, VOLT={}, LAG={}, L={}, FREQ={};'.format(
				self.name, self.VOLT, self.LAG, self.L, self.FREQ)

class DriftDef(Definition):
	def __init__(self, name, L):
		self.name = name
		self.L = L

	def toString(self):
		return '{}: DRIFT, L={};'.format(self.name, self.L)

class SBendDef(Definition):
	def __init__(self, name, L, ANGLE, TILT, K1, K1S, K2, THICK=False, FINT=0, FINTX=0):
		self.name = name
		self.L = L
		self.ANGLE = ANGLE
		self.TILT = TILT
		self.K1 = K1
		self.K1S = K1S
		self.K2 = K2
		self.THICK = THICK
		self.FINT = FINT
		self.FINTX = FINTX

	def toString(self):
		return '{}: SBEND, L={}, ANGLE={}, TILT={}, K1={}, K1S={}, K2={}, THICK={}, FINT={}, FINTX={};'.format(
			self.name, self.L, self.ANGLE, self.TILT, self.K1, self.K1S, self.K2, self.THICK, self.FINT, self.FINTX)

class DipoleEdgeDef(Definition):
	def __init__(self, name, H, E1, FINT, HGAP, TILT):
		self.name = name
		self.H = H
		self.E1 = E1
		self.FINT = FINT
		self.HGAP = HGAP
		self.TILT = TILT


	def toString(self):
		return '{}: DIPEDGE, H={}, E1={}, FINT={}, HGAP={}, TILT={};'.format(self.name, self.H,
			self.E1, self.FINT, self.HGAP, self.TILT)

class NLLensDef(Definition):
	def __init__(self, name, KNLL, CNLL):
		self.name = name
		self.KNLL = KNLL
		self.CNLL = CNLL

	def toString(self):
		return '{}: NLLENS, KNLL={}, CNLL={};'.format(self.name, self.KNLL, self.CNLL)

class MatrixDef(Definition):
	def __init__(self, name, L, KICK, RM, TM):

		self.name = name
		self.L = L
		self.KICK = ''
		self.RM = ''
		self.TM = ''

		for i in range(6):
			self.KICK += 'KICK{}={},'.format(i+1, KICK[i])
			for k in range(6):
				self.RM += 'RM{}{}={},'.format(i+1, k+1, RM[i,k])
		self.KICK = self.KICK[:-1]
		self.RM = self.RM[:-1]
		self.TM = 0 #self.TM[:,-1]

	# TM is not implemented since it cannot be exported to SixTrack
	def toString(self):
		return '{}: MATRIX, L={}, {}, {};'.format(self.name, self.L,
			self.KICK, self.RM)

class MarkerDef(Definition):
	def __init__(self, name):
		self.name = name

	def toString(self):
		return '{}: MARKER;'.format(self.name)

class Lattice:
	def __init__(self, circum):
		self.circum = circum 	 # Circumference of the accelerator
		self.definitions=[]	 	 # List of Definition objects
		self.elements=[]		 # List of tuple(defname, pos)

	# Add a defintion of a "pure" multipole to the lattice
	def addMultipoleDef(self, name, order, KN, KS=0, TILT=0, THICK=False):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = MultipoleDef(name, order, KN, KS, TILT, THICK)
		self.definitions.append(d)

	#a Add a definition of a solenoid to the lattice
	def addSolenoidDef(self, name, L, KS, KSI):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = SolenoidDef(name, L, KS, KSI)
		self.definitions.append(d)

	# Add a definition of a "pure" RF Multipole to the lattice
	def addRFMultipoleDef(self, name, order, VOLT, LAG, HARMON,
		FREQ, TILT, KN, KS, PN, PS):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = RFMultipoleDef(name, order, VOLT, LAG, HARMON, FREQ, TILT,
			KN, KS, PN, PS)
		self.definitions.append(d)

	# Add a definition of a Crab Cavity to the lattice
	def addCrabCavityDef(self, name, L, TILT, VOLT, LAG, HARMON, FREQ, RV1,
			RV2, RV3, RV4, LAGF, RPH1, RPH2):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = CrabCavityDef(name, L, TILT, VOLT, LAG, HARMON, FREQ, RV1,
			RV2, RV3, RV4, LAGF, RPH1, RPH2)
		self.definitions.append(d)

	# Add a definition of a Beam-Beam element to the lattice
	def addBeamBeamDef(self, name, CHARGE, XMA, YMA, SIGX, SIGY, WIDTH, BBSHAPE, BBDIR):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = BeamBeamDef(name, CHARGE, XMA, YMA, SIGX, SIGY, WIDTH, BBSHAPE, BBDIR)
		self.definitions.append(d)

	# Add a definition of a drift space to the lattice
	def addDriftDef(self, name, L):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = DriftDef(name, L)
		self.definitions.append(d)

	# Add a defintion of a bending magnet (type S) to the lattice
	def addSBendDef(self, name, L, ANGLE, TILT, K1, K1S, K2, THICK=False, FINT=0, FINTX=0):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = SBendDef(name, L, ANGLE, TILT, K1, K1S, K2, THICK, FINT, FINTX)
		self.definitions.append(d)

	# Add a defintion of a dipole edge to the lattice
	def addDipoleEdgeDef(self, name, H, E1, FINT, HGAP, TILT):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = DipoleEdgeDef(name, H, E1, FINT, HGAP, TILT)
		self.definitions.append(d)

	# Add a defintion of a Kicker to the lattice
	def addKickerDef(self, name, HKICK, VKICK, L, TILT):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = KickerDef(name, HKICK, VKICK, L, TILT)
		self.definitions.append(d)	
	
	# Add a definition of an RF Cavity to the lattice
	def addRFCavityDef(self, name, VOLT, LAG, L, HARMON, FREQ):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = RFCavityDef(name, VOLT, LAG, L, HARMON, FREQ)
		self.definitions.append(d)

	# Add a definition of a nonlinear lens w/ elliptic potential
	# NOT IMPLEMENTED IN SIXTRACK
	def addNLLensDef(self, name, KNLL, CNLL):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = NLLensDef(name, KNLL, CNLL)
		self.definitions.append(d)

	# Add a definition of an Arbitrary Matrix Element
	def addMatrixDef(self, name, L, KICK, RM, TM=0):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = MatrixDef(name, L, KICK, RM, TM)
		self.definitions.append(d)

	# Add a definition of a Marker
	def addMarkerDef(self, name):
		self.definitions = [d for d in self.definitions if d.name != name]
		d = MarkerDef(name)
		self.definitions.append(d)

	# Add element previously defined by identifier name at position pos
	def addElement(self, name, pos):
		if name in [d.name for d in self.definitions]:
			if (name, pos) in self.elements:
				self.elements.remove((name, pos))
			self.elements.append((name, pos))

	# Return the lattice as a tuple of two strings and the circumference
	def getLatticeDefinition(self, sort=1):
		# Construct string for definitions
		s_defs = ''
		for d in self.definitions:
			s_defs += d.toString() + '\n'

		# Sort elements w.r.t. position
		if sort == 1:
			self.elements.sort(key=lambda tup : tup[1])

		# Construct string for elements
		s_elems = ''
		for n, pos in self.elements:
			s_elems += '{}: {}, at = {};\n'.format(n, n, pos)

		return s_defs[:-1], s_elems[:-1], self.circum

# A basic lattice used as reference for some simple testing
def getReferenceLattice():
	lattice = Lattice(20)
	lattice.addMultipoleDef(name='qf', order=2, KN=0.11755705)
	lattice.addMultipoleDef(name='qd', order=2, KN=-0.11755705)
	lattice.addRFCavityDef(name='cav', VOLT=100, LAG=0.0, L=0.0, HARMON=100, FREQ=0)
	lattice.addElement(name='qd', pos=10)
	lattice.addElement(name='cav', pos=10.001)
	lattice.addElement(name='qf', pos=19.9999)
	return lattice
