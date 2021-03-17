import math
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator
from scipy.optimize import curve_fit
import numpy as np

iFilename = 'cry_interaction.dat'

# Reading file -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with open(iFilename,'r') as readFile:
	for line in readFile:
		partID, turn, collimator, prev, proc, kickx, kicky, Ein, Eout, xpin, ypin, cryangle, xin, yin = np.genfromtxt(readFile, unpack=True)

# Building histograms ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

binx = 200
biny = 200
binang = 200
pmin = -100e-6
pmax = 200e-6
kickmin = -50e-6
kickmax = 100e-6

f1 = plt.figure(1)
plt.hist2d(xpin, kickx, bins=[binang, binang], range=[(pmin,pmax), (kickmin,kickmax)], norm=LogNorm())
plt.title('Deflection along x')
plt.xlabel('x\' in [rad]')
plt.ylabel('kick x [rad]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.colorbar()

f1.show()
f1.show()
f1.set_size_inches((9, 6), forward=False)
f1.tight_layout(rect=[0, 0.03, 1, 0.95])
f1.savefig('deflectionx.png')

f2 = plt.figure(2)
plt.hist2d(ypin, kicky, bins=[binang, binang], range=[(pmin,pmax), (kickmin,kickmax)], norm=LogNorm())
plt.title('Deflection along ya')
plt.xlabel('y\' in [rad]')
plt.ylabel('kick y [rad]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.colorbar()

f2.show()
f2.set_size_inches((9, 6), forward=False)
f2.tight_layout(rect=[0, 0.03, 1, 0.95])
f2.savefig('deflectiony.png')

input()
plt.close('all')

