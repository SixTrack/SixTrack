import compareSix2Mad as csm
import latticeConstructor as lc
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
import pandas as pd
import os

# ------ Documentation on output format for tests ------
# n : an ndarray of shape (nbr_of_test_cases, points_per_test_case) containing parameter values
#     For the ith test case, the values of the parameter being iteratively altered
#     is given by n(i-1,:).
# vals : an ndarray of shape (6, nbr_of_test_cases, 3, points_per_test_case) containing test output
#        The specifics of each dimension:
#        - The 1st dimension of vals is for the six different MAD-X tracking variables in order (X,PX,Y,PY,T,PT)
#        - The 2nd dimension of vals is for each different test case
#        - The 3rd dimension of vals is for the MAD-X output, the SixTrack output in MAD-X format and
#          the difference between them respectively
#        - The 4th dimension of vals is for each parameter value
#
# Examples:
#           - vals(1,2,0,5) corresponds to the PX value in MAD-X for the third test case in which the parameter
#             value is n(2,5)
#           - vals(4,0,2,7) corresponds to the difference in T value between MAD-X and MAD-X-transformed SixTrack for
#             the first test case in which the parameter value is n(0,7)
#           - vals(2,1,1,1) corresponds to the Y value in MAD-X-transformed SixTrack for the second test case in which
#             the parameter value is n(1,1)
#

nbr_turns = 1
energy_GeV = 5
init_coords = (0.001, 0.002, 0.003, 0.004, 0.005, 0.01)
rl = lc.getReferenceLattice()
path = os.getcwd()

def loadValues(vals, i, j, diff, mad, six2mad):
    for col_idx in range(6):
        vals[col_idx, i, 0, j] = mad.iloc[0, 2+col_idx]
        vals[col_idx, i, 1, j] = six2mad.iloc[0, col_idx]
        vals[col_idx, i, 2, j] = diff.iloc[0, col_idx]

def testMultipole(orde):
    rl = lc.getReferenceLattice()

    vals = np.zeros((6,4,3,100))
    n = np.zeros((4,100))
    
    n[0,:] = [0.001 * x for x in range(1,101)]
    for j in range(100):
        rl.addMultipoleDef(name='trkstrt', order=orde, KN= (j+1) * 0.0008, KS=0, TILT=0, THICK=False)
        rl.addElement(name='trkstrt', pos=0.0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 0, j, diff, mad, six2mad)

    n[1,:] = [0.001 * x for x in range(1,101)]
    for j in range(100):
        rl.addMultipoleDef(name='trkstrt', order=orde, KN= 0, KS=(j+1) * 0.0008, TILT=0, THICK=False)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 1, j, diff, mad, six2mad)

    n[2,:] = [np.pi * 0.02 * x for x in range(100)]
    for j in range(100):
        rl.addMultipoleDef(name='trkstrt', order=orde, KN= 0.008, KS=0, TILT= str(0.02*j)+'*PI',
            THICK=False)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 2, j, diff, mad, six2mad)

    n[3,:] = [np.pi * 0.02 * x for x in range(100)]
    for j in range(100):
        rl.addMultipoleDef(name='trkstrt', order=orde, KN= 0, KS=0.008, TILT=str(0.02*j)+'*pi',
            THICK=False)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)
        
        loadValues(vals, 3, j, diff, mad, six2mad)
        
    return n, vals

def testMatrix():
    rl = lc.getReferenceLattice()
    kick = np.zeros(6)
    rm = np.eye(6)
    
    vals = np.zeros((6,6,3,100))
    n = np.zeros((6,100))

    for k in range(6):
        n[k,:] = [0.0001 * (x+1) for x in range(100)]

    for i in range(6):
        for j in range(100):
            kick = np.zeros(6)
            kick[i] = (j+1) * 0.0001

            rl.addMatrixDef(name='trkstrt', L=0.0, KICK=kick, RM=rm)
            rl.addElement(name='trkstrt', pos=0)
            rl_s = rl.getLatticeDefinition()
            
            diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
                rl_s, norm='', verbose=1, all_files=1)

            loadValues(vals, i, j, diff, mad, six2mad)
    return n, vals

def driftTest():
    rl = lc.getReferenceLattice()
    
    vals = np.zeros((6,1,3,100))
    n = np.zeros((1,100))

    n[0,:] = [0.1 * x for x in range(100)]
    for j in range(100):
        rl = lc.getReferenceLattice()
        rl.addMarkerDef('IP1')
        rl.addElement(name='IP1', pos= j * 0.1)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('IP1', 'ip1'))

        loadValues(vals, 0, j, diff, mad, six2mad)

    return n, vals

# NOTE: MAD-X appends a 'd/q/s/o' and  an 's' for KS!=0 
# to the element name of an exported RFMultipole, hence the
# need for the symb tuple
# NOTE: RF Sextupole crashes this test (as in does not run to completion). The element itself is verified
# but in the MAD-X to SixTrack conversion it produces a misalignment regardless of if it is misaligned,
# and it is wrongly named (i.e. not appended w/ a 's') which ruins the execution of the test below. Not
# copying fc.8 to fort.8 unless it is actually misaligned allows the test to be run contrary to this.
def testRFMultipole(orde):
    symb = ('', 'd', 'q', 's', 'o')
    rl = lc.getReferenceLattice()

    vals = np.zeros((6,7,3,100))
    n = np.zeros((7,100))
    
    n[0,:] = [0.0002 * x for x in range(1,101)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG=0.25, HARMON=0,
            FREQ=200, TILT=0, KN=0.0002 * (j+1), KS=0.0, PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0.0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]))

        loadValues(vals, 0, j, diff, mad, six2mad)

    n[1,:] = [0.0002 * x for x in range(1,101)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG=0.25, HARMON=0,
            FREQ=200, TILT=0, KN=0.0, KS=0.0002 * (j+1), PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0.0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]+'s'))

        loadValues(vals, 1, j, diff, mad, six2mad)

    n[2,:] = [2*(x+1) for x in range(100)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG=0.25, HARMON=0,
            FREQ=2*(j+1), TILT=0, KN=0.02, KS=0.0, PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]))
        
        loadValues(vals, 2, j, diff, mad, six2mad)

    n[3,:] = [x * 0.01 for x in range(100)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG=0.25, HARMON=0,
            FREQ=200, TILT=0, KN=0.02, KS=0.0, PN=0.01 * j, PS=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]))

        loadValues(vals, 3, j, diff, mad, six2mad)

    n[4,:] = [np.pi * x / 50 for x in range(100)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG= 0.25,
            HARMON=0, FREQ=200, TILT=np.pi * j /50, KN=0.02, KS=0.0, PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]))

        loadValues(vals, 4, j, diff, mad, six2mad)


    n[5,:] = [np.pi * x / 50 for x in range(100)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0, LAG= 0.25,
            HARMON=0, FREQ=200, TILT=np.pi * j /50, KN=0.0, KS=0.02, PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]+'s'))

        loadValues(vals, 5, j, diff, mad, six2mad)

    n[6,:] = [0.02*x for x in range(100)]
    for j in range(100):
        rl.addRFMultipoleDef(name='trkstrt', order=orde, VOLT=0.02*j, LAG= 0.25,
            HARMON=0, FREQ=200, TILT=0, KN=0.02, KS=0.0, PN=0.0, PS=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()
        
        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1, trk_element=('trkstrt', 'trkstrt'+symb[orde]))

        loadValues(vals, 6, j, diff, mad, six2mad)

    return n, vals

def testRFCavity():
    rl = lc.getReferenceLattice()

    vals = np.zeros((6,3,3,100))
    n = np.zeros((3,100))

    n[0,:] = range(1,101)
    for j in range(100):

        rl.addRFCavityDef(name='trkstrt', VOLT=j+1, LAG=0.25, L=0, HARMON=100, FREQ=0)
        rl.addElement(name='trkstrt', pos=0.0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 0, j, diff, mad, six2mad)

    n[1,:] = [0.01 * x for x in range(100)]
    for j in range(100):
        rl.addRFCavityDef(name='trkstrt', VOLT=100, LAG= 0.01*j, L=0, HARMON=100, FREQ=0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 1, j, diff, mad, six2mad)

    n[2,:] = range(1,101)
    for j in range(100):
        rl.addRFCavityDef(name='trkstrt', VOLT=100, LAG=0.25, L=0, HARMON=j+1, FREQ=0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 2, j, diff, mad, six2mad)

    return n, vals

def testSolenoid():
    rl = lc.getReferenceLattice()

    vals = np.zeros((6,3,3,100))
    n = np.zeros((3,100))
    
    n[0,:] = [0.1 * x for x in range(100)]
    for j in range(100):

        rl.addSolenoidDef(name='trkstrt', KS = j * 0.1, KSI = 0.000001, L=0.0)
        rl.addElement(name='trkstrt', pos=0.0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 0, j, diff, mad, six2mad)

    n[1,:] = [0.01 * x for x in range(100)]
    for j in range(100):
        rl.addSolenoidDef(name='trkstrt', KS = 0.000001, KSI = j * 0.01, L=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 1, j, diff, mad, six2mad)

    n[2,:] = [0.005 * x for x in range(100)]
    for j in range(100):
        rl.addSolenoidDef(name='trkstrt', KS = j * 0.005, KSI = j * 0.005, L=0.0)
        rl.addElement(name='trkstrt', pos=0)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 2, j, diff, mad, six2mad)
    return n, vals

def testKicker():
    rl = lc.getReferenceLattice()

    vals = np.zeros((6,3,3,100))
    n = np.zeros((3,100))

    n[0,:] = [0.0003 * (x+1) for x in range(100)]
    for j in range(100):
        rl.addKickerDef(name='trkstrt', HKICK= 0.0003 * (j+1), VKICK=0, L=0, TILT=0)
        rl.addElement(name='trkstrt', pos=0.01)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 0, j, diff, mad, six2mad)


    n[1,:] = [0.0003 * (x+1) for x in range(100)]
    for j in range(100):
        rl.addKickerDef(name='trkstrt', HKICK= 0, VKICK=0.0003 * (j+1), L=0, TILT=0)
        rl.addElement(name='trkstrt', pos=0.01)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 1, j, diff, mad, six2mad)

    n[2,:] = [np.pi * x / 50 for x in range(100)]
    for j in range(100):
        rl.addKickerDef(name='trkstrt', HKICK= 0.003, VKICK=0, L=0, TILT=np.pi * j / 50)
        rl.addElement(name='trkstrt', pos=0.01)
        rl_s = rl.getLatticeDefinition()

        diff, mad, six2mad = csm.compare('element_test', nbr_turns, energy_GeV, init_coords, 0,
            rl_s, norm='', verbose=1, all_files=1)

        loadValues(vals, 2, j, diff, mad, six2mad)

    return n, vals

if len(argv) == 1:
    print('Enter an option!')
    exit()
else:
    if argv[1] == 'rfcav':
        n, vals = testRFCavity()
    elif argv[1] == 'multipole':
        if len(argv) == 2:
            n, vals = testMultipole(1)
        else:
            n, vals = testMultipole(int(argv[2]))
    elif argv[1] == 'kicker':
        n, vals = testKicker()
    elif argv[1] == 'solenoid':
        n, vals = testSolenoid()
    elif argv[1] == 'rfmultipole':
        if len(argv) == 2:
            n, vals = testRFMultipole(1)
        else:
            n, vals = testRFMultipole(int(argv[2]))
    elif argv[1] == 'matrix':
        n, vals = testMatrix()
    elif argv[1] == 'drift':
        n, vals = driftTest()

print(n)
print(vals)