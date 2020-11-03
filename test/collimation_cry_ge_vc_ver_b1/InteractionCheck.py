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

# Calculations -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

countAM = 0
countCH = 0
countVR = 0
countVC = 0
total = len(partID)

for i in range(total):
	if proc[i] == 1:
		countAM = countAM+1
	if proc[i] == 2:
		countVR = countVR+1
	if proc[i] == 3:
		countCH = countCH+1
	if proc[i] == 4:
		countVC = countVC+1

print('--- Interactions at Cry ---')
print('Amorphous: ', countAM, "{0:.2%}".format(countAM/float(total)))
print('Volume reflection: ', countVR, "{0:.2%}".format(countVR/float(total)))
print('Volume capture: ', countVC, "{0:.2%}".format(countVC/float(total)))
print('Channeling: ', countCH, "{0:.2%}".format(countCH/float(total)))

