import math
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator
from scipy.optimize import curve_fit
import numpy as np

iFilename = 'cry_exit.dat'

# Reading file -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with open(iFilename,'r') as readFile:
	for line in readFile:
		partID, turn, collimator, mat, ishit, isabs, x, xp, y, yp, p = np.genfromtxt(readFile, unpack=True)

# Building histograms ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

binx = 200
biny = 200
binang = 200
xmin = -3e-3
xmax = 3e-3
xpmin = -3e-5
xpmax = 9e-5
ymin = -8e-4
ymax = 8e-4
ypmin = -1.5e-5
ypmax = 1.5e-5

f0 = plt.figure(0)
plt.suptitle('Beam profile @ crystal exit')

plt.subplot(2,2,1)
# plt.hist(x, bins=binx, range=(xmin,xmax))
plt.hist(x, bins=binx)
plt.title('x profile')
plt.xlabel('x [m]')
plt.ylabel('# impacts')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale('log')

plt.subplot(2,2,2)
# plt.hist(y, bins=biny, range=(ymin,ymax))
plt.hist(y, bins=biny)
plt.title('y profile')
plt.xlabel('y [m]')
plt.ylabel('# impacts')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale('log')

plt.subplot(2,2,3)
# plt.hist(xp, bins=binang, range=(xpmin,xpmax))
plt.hist(xp, bins=binang)
plt.title('x\' profile')
plt.xlabel('x\' [rad]')
plt.ylabel('# impacts')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale('log')

plt.subplot(2,2,4)
# plt.hist(yp, bins=binang, range=(ypmin,ypmax))
plt.hist(yp, bins=binang)
plt.title('y\' profile')
plt.xlabel('y\' [rad]')
plt.ylabel('# impacts')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.yscale('log')

mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())

f0.show()
f0.set_size_inches((13, 9), forward=False)
f0.tight_layout(rect=[0, 0.03, 1, 0.95])
f0.savefig('profiles1d_exit.png')

f1 = plt.figure(1)
# plt.hist2d(x, y, bins=[binx, biny], range=((xmin,xmax),(ymin,ymax)), norm=LogNorm())
plt.hist2d(x, y, bins=[binx, biny], norm=LogNorm())
plt.title('2d profile @ crystal exit')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
# plt.axvline(x=collgap*1.0e3)
# plt.axvline(x=mincoll*1.0e3)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.colorbar(ticks = LogLocator(subs=range(10)))

f1.show()
f1.show()
f1.set_size_inches((9, 6), forward=False)
f1.tight_layout(rect=[0, 0.03, 1, 0.95])
f1.savefig('profile2d_exit.png')

f2 = plt.figure(2)
# plt.hist2d(y, yp, bins=[biny, binang], range=((xmin,xmax),(xpmin,xpmax)), norm=LogNorm())
plt.hist2d(y, yp, bins=[biny, binang], norm=LogNorm())
plt.title('y-y\' phase space @ crystal exit')
plt.xlabel('y [m]')
plt.ylabel('y\' [rad]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.colorbar(ticks = LogLocator(subs=range(10)))

f2.show()
f2.set_size_inches((9, 6), forward=False)
f2.tight_layout(rect=[0, 0.03, 1, 0.95])
f2.savefig('phasespacey_exit.png')

raw_input()
plt.close('all')

