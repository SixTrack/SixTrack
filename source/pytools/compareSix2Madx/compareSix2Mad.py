import numpy as np
import pandas as pd
import os,shutil
speedOfLight = 299792458;
massProton=0.938272081

def energyToMomentum(energy):
    print energy, np.sqrt(energy**2-massProton**2)
    return np.sqrt(energy**2-massProton**2)

def momentumToEnergy(momentum):
    return np.sqrt(momentum**2+massProton**2)

def getRelBeta(energy):
    return (energyToMomentum(energy)/energy)

def getSigmaZfromT(energy, t):
    return getRelBeta(energy)*t

def getPxFromAngle(angle, de):
    return angle/(1+de)

def getPTfromDeltaP(energy_GeV,delta_p):
    p0 = energyToMomentum(energy_GeV)
    pc = (1+delta_p)*p0
    energy_withDpGeV= momentumToEnergy(pc)
    return (energy_withDpGeV-energy_GeV)/p0
def getDeltaSFromT(energy0, time, pt):
    return 0
def getBetaOverBeta0(energy, deltap):
    return 0

def createMultiAndPhase(order,strenght, phase):
    knl = ''
    ph = ''
    for i in range(0, order-1):
        knl = knl +'0, '
        ph = ph +'0, '
    
    knl = knl + str(strenght)
    ph = ph + str(phase)
    return [knl, ph]

if(0 is 1):
    path = '/afs/cern.ch/work/t/tpersson/simpleTrack/onlyonePASS/'
    printToFile=1
    nameFile = 'compare.csv'
    madxPath = '/home/tobias/codes/MAD-X/madx64  '
    sixtrackPath = './newTry '
    #Parameters for the simulation
    energy_GeV = 2
    deltap = 0.1
    x_m=0.002
    y_m=0.003
    delta_p=0.1
    time_mad =0.1

    #Convert the parameters to corresponding output
    x_mm = x_m*1000
    y_mm = y_m*1000
    energy_MeV=energy_GeV*1000
    pt_mad = getPTfromDeltaP(energy_GeV,delta_p) 
    sigma_z = getSigmaZfromT(energy_GeV, time_mad)*1000

    #Set up fort.3 
    intext=file(path+ 'fort.macro.3','r').read()
    file(path+'fort.3','w').write(intext % locals())


    multi = range(1,5)
    arrPandas= []
    diffSix=  pd.DataFrame()
    diffSix['xp']=0

    #RF phases
    phases = [0.001, 0.12, 0.25]
    k_value = 0.001
    for mult in multi:
        for phase in phases:
            
            [knl_strength, knl_phase] = createMultiAndPhase(mult, k_value, phase)
            
            intext=file(path+ 'job.macro.madx','r').read()
            file(path+'job.madx','w').write(intext % locals())
            os.chdir(path)
            
            os.system(madxPath + 'job.madx')
            os.system('cp fc.2 fort.2')
            os.system(sixtrackPath)  
            sixPath =  path + 'crab_dump.txt'
            madPath =  path + 'qfstart.obs0002.p0001'
            six = pd.read_table(sixPath, delim_whitespace =True, names=['ID', 'turn', 's', 'x', 'xp', 'y', 'yp', 'z','de', 'ktrack'], skiprows=2)
            mad = pd.read_table(madPath, delim_whitespace =True, names=['ID', 'turn', 'x', 'xp', 'y', 'yp', 't','pt', 's','E'], skiprows=8)
            
            six_p1 = six.loc[six['ID'] == 1]
            six_p1.index=range(0,len(six_p1))
            
            energy = mad['E'].values[0]
            diffSix=  pd.DataFrame() 
            p0 = energyToMomentum(energy)
            energy_real = energy*(1+six_p1['de'].values[0])
            pc = np.sqrt(energy_real**2-massProton**2)
            dp = (pc-p0)/p0
            six_p1['de'].values
            diffSix['xp'] =  0.001*six_p1['xp'].values*(1+dp)-mad['xp'].values
            diffSix['yp'] =  0.001*six_p1['yp'].values*(1+dp)-mad['yp'].values
            diffSix['x'] =  0.001*six_p1['x'].values-mad['x'].values
            diffSix['y'] =  0.001*six_p1['y'].values-mad['y'].values
            diffSix['pt'] = six_p1['de'].values*energy/energyToMomentum(energy)-mad['pt']
            diffSix['z'] = 0.001*six_p1['z']-getSigmaZfromT(energy, mad['t'].values)
            diffSix['multi'] = mult
            diffSix['phase'] = phase
            arrPandas.append(diffSix)


    if(printToFile):
        if(os.path.isfile(nameFile)):
            os.remove(nameFile)
        for diff in arrPandas:
            if not os.path.isfile(nameFile):
                diff.to_csv(nameFile,header ='column_names')
            else:
                diff.to_csv(nameFile,mode = 'a',header=False)
        tmp=pd.read_csv(path+nameFile)
        print tmp
