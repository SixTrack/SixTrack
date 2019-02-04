import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os,shutil,time

speedOfLight = 299792458;
massProton=0.938272081

# These paths need to be manually set
#   madxPath = absolute path to the madx64 executable
#   sixtrackPath = absolute path to the SixTrack executable

madxPath = '/home/tobias/codes/MAD-X/madx64'
sixtrackPath = '/home/tobias/codes/SixTrackTobias/build/SixTrack_cmakesix_BUILD_TESTING_defaultcompiler_defaultbuildtype/SixTrack_50002_crlibm_rn_Linux_gfortran_static_x86_64_64bit_double'

# E -> P
def energyToMomentum(E):
    return np.sqrt(E**2-massProton**2)

# P -> E
def momentumToEnergy(p):
    return np.sqrt(p**2+massProton**2)

# E -> \beta
def getRelBeta(E):
    return (energyToMomentum(E)/E)

# E0, dE/E0 -> dp/p0
def getDeltaP(E0, dE_frac):
	E = (1+dE_frac) * E0
	p = energyToMomentum(E)
	p0 = energyToMomentum(E0)
	return (p - p0) / p0

# E0, PT -> dp/p0
def getDeltaP2(E0, PT):
    p0 = energyToMomentum(E0)
    E = PT * p0 + E0
    p = energyToMomentum(E)
    return (p-p0)/p0

# E0, T -> \sigma_z
def getSigmaZ(E0, t):
    return getRelBeta(E0) * t * 1000

# E0, \sigma_z -> T
def getT(E0, z):
    return z * 0.001 / getRelBeta(E0)

# E0, dE/E0 -> PT
def getPT(E0, dE_frac):
    return dE_frac * E0 / energyToMomentum(E0)

# E0, PT -> dE/E0
def getDE_frac(E0, PT):
    p0 = energyToMomentum(E0)
    E = PT * p0 + E0
    return (E-E0)/E0

# E0, dE/E0, PX or PY -> x' or y' 
def getAngle(E0, dE_frac, PX):
    return 1000 * PX / (1 + getDeltap(E0, dE_frac))

# dp/p0, PX or PY -> x' or y'
def getAngle2(deltap, PX):
    return 1000 * PX / (1 + deltap)

# x' or y' -> PX or PY
def xpToPX(E0, dE_frac, xp):
    return (1 + getDeltaP(E0, dE_frac)) * (xp / 1000)

# Transforms initial coordinates in MAD-X to initial coordinates for SixTrack
def madToSixCoords(mad_coords, E0_GeV):
    x_m, PX, y_m, PY, time_mad, pt_mad = mad_coords
    x_mm = x_m*1000
    y_mm = y_m*1000
    sigma_z = getSigmaZ(E0_GeV, time_mad)
    delta_p = getDeltaP2(E0_GeV, pt_mad)
    xp = getAngle2(delta_p, PX)
    yp = getAngle2(delta_p, PY)
    E0_MeV = E0_GeV * 1000
    return x_mm, xp, y_mm, yp, sigma_z, delta_p, E0_MeV

# Given a SixTrack tracking output file transforms it into MAD-X 'standard' tracking output
def sixToMad(six, E0):
    dE_frac = six['de'][0::2].values
    z = six['z'][0::2].values
    xp = six['xp'][0::2].values
    yp = six['yp'][0::2].values

    six2mad = pd.DataFrame()
    six2mad['X'] = six['x'][0::2].values / 1000.0
    six2mad['PX'] = xpToPX(E0, dE_frac, xp)
    six2mad['Y'] = six['y'][0::2].values / 1000.0
    six2mad['PY'] = xpToPX(E0, dE_frac, yp)
    six2mad['T'] = getT(E0, z)
    six2mad['PT'] = getPT(E0, dE_frac)
    return six2mad

# Retrieves the difference between two tables containing MAD-X 'standard' variables
def getTableDifference(mad1, mad2, norm):
    diff =  pd.DataFrame()
    if norm == '':
        diff['dX'] =  mad1['X']-mad2['X']
        diff['dPX'] =  mad1['PX']-mad2['PX']
        diff['dY'] =  mad1['Y']-mad2['Y']
        diff['dPY'] =  mad1['PY']-mad2['PY']
        diff['dT'] = mad1['T']-mad2['T']
        diff['dPT'] = mad1['PT']-mad2['PT']
    elif norm == 'abs':
        diff['dX'] =  abs(mad1['X']-mad2['X'])
        diff['dPX'] =  abs(mad1['PX']-mad2['PX'])
        diff['dY'] =  abs(mad1['Y']-mad2['Y'])
        diff['dPY'] =  abs(mad1['PY']-mad2['PY'])
        diff['dT'] = abs(mad1['T']-mad2['T'])
        diff['dPT'] = abs(mad1['PT']-mad2['PT'])
    return diff

# Loads a table with filepath 'path'.
# If names is not given, pandas will try to automatically parse column names from the input
def loadTable(path, isMadX=True, names=[]):
    if isMadX:
        with open(path,'r') as f:
            s = f.read()
        path_temp = path + '_temp'
        with open(path_temp,'w') as f:
            s = s.split('*',1)[1]
            if len(names) == 0:
                s = '\n'.join(s.split('\n', 2)[0:3:2])
            else:
                s = s.split('\n', 2)[2]
            f.write(s)
        if len(names) == 0:
            table = pd.read_table(path_temp, delim_whitespace=True)
        else:
            table = pd.read_table(path_temp, delim_whitespace=True, names=names)
        os.remove(path_temp)
    else:
        if len(names) == 0:
            table = pd.read_table(path, delim_whitespace=True, skiprows=2)
        else:
            table = pd.read_table(path, delim_whitespace=True, names=names, skiprows=2)

    return table

# Retrieves tables from a list of output tables in /work
def getRawTables(file_list):
    path = os.getcwd()
    table_list = []
    for filename, isMadX, names in file_list:
        filepath = path + '/work/' + filename
        table_list.append(loadTable(filepath, isMadX, names))
    return table_list

# Description:
# Takes two (absolute or relative) paths to the tracking output of a SixTrack and MAD-X run, compares
# the values in MAD-X format and returns the result
#
# Parameters:
#   mad_output_path : Absolute path to MAD-X output file
#   six_output_path : Absolute path to SixTrack output file
#   E0_GeV : Energy of the reference particle in GeV
#   all_files : 1 means that the a dataframe of the difference, the converted SixTrack output and the
#               MAD-X output are returned, otherwise only the difference is returned
#   printToFile : 1 means a file of called [filename] gets printed in output/
#   filename : name of the output diff file
#   norm : defines how the difference is to be computed. Currently the following
#          values are legal
#           - '' : take the signed difference
#           - 'abs' : take the absolute value of the difference
#
def compareOutput(mad_output_path, six_output_path, E0_GeV, all_files=0, printToFile=0,
    filename='diff.csv', norm=''):

    # Get current working directory
    path = os.getcwd()

    # Read the output
    mad = loadTable(mad_output_path, names=['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T','PT', 'S','E'],
    	isMadX=True)
    six = loadTable(six_output_path, names=['ID', 'turn', 's', 'x', 'xp', 'y', 'yp', 'z','de', 'ktrack'],
    	isMadX=False)

    # Transform result in SixTrack to MAD-X standard
    six2mad = sixToMad(six, E0_GeV)

    # Construct table of differences in tracking
    diff = getTableDifference(six2mad, mad, norm)

    # Print difference to file in compareSix2Mad/output directory
    if printToFile:
        diff.to_csv('work/' + filename, header = 'column_names')
    if all_files == 1:
        return diff, mad, six2mad
    else:
        return diff


# Run MAD-X and SixTrack for identical lattice + initial conditions inside the /work directory
# CAUTION: WILL DELETE ALL CONTENT IN THE CURRENT /work DIRECTORY IN THE PROCESS
def runSimulations(test_name, mad_coords, nbr_turns, E0_GeV, lattice=(-1,0,0),
        trk_element=('trkstrt', 'trkstrt'), verbose=0):
    
    # Set path to current working directory (should AND NEEDS to be comparSix2Mad directory)
    path = os.getcwd()

    # NOTE: it is possible that the time resolution is insufficient granted the simulations and I/O run too fast
    start_time = time.time()

    # Load the lattice used for the run, if present
    if lattice[0] == -1:
        if verbose == 1:
            print('\nNo lattice loaded, using what is predefined in *.madx\n')
    else:
        s_defs, s_elems, circum = lattice

    # Load name of element from which tracking is to be executed
    mad_track_element, six_track_element = trk_element

    # Load MAD-X paramters
    x_m, PX, y_m, PY, time_mad, pt_mad = mad_coords

    # Clean or create the /work directory
    if os.path.isdir(path + '/work'):
        shutil.rmtree(path + '/work')
    os.mkdir(path + '/work')

    #Set up job.madx
    with open(path + '/scenarios/' + test_name + '.macro.madx','r') as f_macro:
        data = f_macro.read()

        with open(path + '/work/' +'job.madx','w') as f:
            f.write(data % locals())

    # Change to /work directory
    os.chdir(path + '/work/')

    # Run MAD-X on job.madx
    if verbose == 1:
        os.system(madxPath + ' job.madx')
    else:
        os.system(madxPath + ' job.madx > /dev/null')

    #Convert the MAD-X parameters to SixTrack parameters
    x_mm, xp, y_mm, yp, sigma_z, delta_p, E0_MeV = madToSixCoords(mad_coords, E0_GeV)

    # Copy exported structure to fort.2
    shutil.copyfile('fc.2', 'fort.2')

    # Read fc.3 and fc.3.aux for input to fort.3 if applicable
    with open('fc.3','r') as f:
        trom_s = f.read()
    with open('fc.3.aux','r') as f:
        data = f.read()
        fc3_aux_s = '\n'.join(data.split('\n')[:7])

    #Set up fort.3
    with open(path + '/scenarios/' + test_name + '.macro.3','r') as f_macro:
        data = f_macro.read()
        with open('fort.3','w') as f:
            f.write(data % locals())

    #Set up fort.13 (if relevant for test)
    macropath = path + '/scenarios/' + test_name + '.macro.13'
    if os.path.isfile(macropath):
        E0_MeV_mod = E0_MeV * (1 + getDE_frac(E0_GeV, pt_mad))
        with open(macropath,'r') as f_macro:
            data = f_macro.read()
            with open('fort.13','w') as f:
                f.write(data % locals())

    # Copy misalignments to fort.8 if relevant
    if os.path.isfile('fc.8'):
        shutil.copyfile('fc.8', 'fort.8')

    # Copy fc.16 to fort.16 if relevant
    if os.path.isfile('fc.16'):
        shutil.copyfile('fc.16', 'fort.16')

    # Run SixTrack on fort.3
    if verbose == 1:
        print('\nFinished running MAD-X\n')
        os.system(sixtrackPath)
        print('\nFinished running SixTrack\n')
    else:
        os.system(sixtrackPath + ' > /dev/null')

    # Switch back to previous cwd before finishing
    os.chdir(path)

    # Now all files should be produced in the /work/ directory
    return


# Description: Runs the same accelerator setup on MAD-X and SixTrack and then outputs
#              the difference between them in MAD-X coordinates. Assumes a specific
#              setup elaborated on in the README in [NOT YET PRESENT]
#
#
# Parameters:
#   test_name : name of MAD-X and SixTrack macro files in scenarios/
#   nbr_turns : number of turns in the simulation
#   E0_GeV : Energy of the refernce particle in GeV
#   mad_coords : initial coordinates given in MAD-X variables (6-tuple)
#   printToFile :" 1 means a file 'diff.csv' gets printed in output/
#   lattice : A tuple defining the lattice of form (a,b,c)
#             where     a : a formatted string containing all the element
#                           definitions for the lattice in MAD-X syntax (see latticeConstructor.py
#                           for an easy way of producing this)
#                       b : a formatted stringing containing the element sequence
#                           defining the lattice in MAD-X syntax (see latticeConstructor.py for an
#                           easy way of producing this)
#                       c : a number (float, integer) defining the circumference
#                           of the accelerator lattice
#   trk_element : A 2-tuple containing the element names from which the particles in MAD-X respectively
#                 SixTrack are to be tracked.
#   norm : defines how the difference is to be computed. Currently the following
#          values are legal
#                           ''      : difference is taken with sign
#                           'abs'   : the absolute value of the difference is taken
#   verbose : set to 1 for MAD-X/SixTrack output in terminal
#   all_files : 1 means that a dataframe of the difference, the converted SixTrack output and the
#               MAD-X output are returned, otherwise only the difference is returned
#
def compare(test_name, nbr_turns, E0_GeV, mad_coords, printToFile=0, lattice=(-1,0,0),
		trk_element=('trkstrt','trkstrt'), norm='', verbose=0, all_files=0):
    
    path = os.getcwd()

    # NOTE: it is possible that the time resolution is insufficient granted the simulations and I/O run too fast
    start_time = time.time()

    # Run simulation on both MAD-X and SixTrack
    runSimulations(test_name, mad_coords, nbr_turns, E0_GeV, lattice, trk_element, verbose)

    # Path to output files
    madPath = path + '/work/mad_output.obs0002.p0001'
    sixPath =  path + '/work/six_output.txt'

    # Make sure files were created
    if os.path.isfile(madPath) and os.path.isfile(sixPath):
        # Compare the the two outputs 
        return compareOutput(madPath, sixPath, E0_GeV, printToFile=printToFile,
            norm=norm, all_files=all_files)
    else:
        # Something went wrong..
        print('Output of SixTrack or MAD-X was not produced during execution.')
        print('Check output of the runs (w/ verbose=1) for more information.')
        return -1
