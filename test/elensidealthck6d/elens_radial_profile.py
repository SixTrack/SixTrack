import matplotlib.pyplot as plt
import numpy as np

sig=4.82065 # 1 sigma beam size, hel1 between 4-6 sigma
mm2cm=0.1


def parseFort6(iFileName='fort.6'):
    'read fort.6 file and get file parsing and normalisation'
    parsedFile=[]
    integratedProfile=[]
    rrs=[]
    print 'parsing %s file...'%(iFileName)
    with open(iFileName,'r') as iFile:
        readProf=False
        readIntegr=False
        for line in iFile.readlines():
            if ( line.startswith('ELENS> Radial profile') ):
                readProf=True
                continue
            elif ( line.startswith('ELENS> Normalising radial profile') ):
                readProf=False
                continue
            elif ( line.startswith('ELENS> Integrated radial profile') ):
                readIntegr=True
                continue
            elif ( line.startswith('ELENS> Total current in radial') ):
                readIntegr=False
                continue
            elif ( line.startswith('ELENS> NB: point at ii=0') ):
                continue
            if ( readProf ):
                data=line.strip().split(',')
                rrs.append(float(data[1]))
                parsedFile.append(float(data[2]))
            elif ( readIntegr ):
                data=line.strip().split(',')
                integratedProfile.append(float(data[2]))
    if (len(rrs)>0):
        print ' ...acquired profile of current in elens from file: %i,%i,%i '% \
            (len(rrs),len(parsedFile),len(integratedProfile))
    else:
        print ' ...strange: no profiles...'
    return np.array(rrs), np.array(parsedFile), np.array(integratedProfile)


rrs, parsedFile, integratedProfile = parseFort6()

profIN=np.loadtxt('CHG1b_170523_8p75A_2-4-2kG_500V_75mA_hires_j-vs-r.txt')
f, axarr = plt.subplots(2,2,figsize=(14,14))
#
axarr[0,0].grid()
axarr[0,0].set_xlabel(r'$r$ [mm]')
axarr[0,0].set_ylabel(r'$J$ [A cm$^{-2}$]')
axarr[0,0].xaxis.set_ticks(np.arange(0,90,10))
if (len(rrs)>0):
    axarr[0,0].plot(rrs,parsedFile*100,'rs-',markeredgewidth=0.0,label='6T echo')
axarr[0,0].plot(profIN[:,0],profIN[:,1],'b.-',label='file')
axarr[0,0].legend(loc='best')
#
axarr[0,1].grid()
axarr[0,1].set_xlabel(r'$r$ [$\sigma$]')
axarr[0,1].set_ylabel(r'$J$ [A cm$^{-2}$]')
axarr[0,1].xaxis.set_ticks(np.arange(0,20,2))
if (len(rrs)>0):
    axarr[0,1].plot(rrs/sig,parsedFile*100,'rs-',markeredgewidth=0.0,label='6T echo')
axarr[0,1].plot(profIN[:,0]/sig,profIN[:,1],'b.-',label='file')
axarr[0,1].legend(loc='best')

# cumulative
cumulIN=[0.0]
tot=0.0
oldVal=0.0
oldR=0.0
intRange=np.append(0.0,profIN[:,0])
for ii in range(len(profIN[:,1])):
    if (ii==0):
        rMin=0.0
    else:
        rMin=profIN[ii-1,0]
    tot+=profIN[ii,1]*np.pi*(profIN[ii,0]-rMin)*(profIN[ii,0]+rMin)*mm2cm**2
    cumulIN.append(tot)
    oldVal=profIN[ii,1]
    oldR=profIN[ii,0]
cumulIN=np.array(cumulIN)
axarr[1,0].grid()
axarr[1,0].set_xlabel(r'$r$ [$\sigma$]')
axarr[1,0].set_ylabel(r'$I$ [A]')
axarr[1,0].xaxis.set_ticks(np.arange(0,20,2))
if (len(rrs)>0):
    axarr[1,0].plot(rrs/sig,integratedProfile,'rs-',markeredgewidth=0.0,label='6T echo')
axarr[1,0].plot(intRange/sig,cumulIN,'b.-',label='file')
axarr[1,0].legend(loc='best')
# natural kick
axarr[1,1].grid()
axarr[1,1].set_xlabel(r'$r$ [$\sigma$]')
axarr[1,1].set_ylabel(r'$\theta$ [A/mm]')
axarr[1,1].xaxis.set_ticks(np.arange(0,20,2))
# rrs[0]=0.0 but it is normal...
if (len(rrs)>0):
    axarr[1,1].plot(rrs/sig,integratedProfile/rrs,'rs-',markeredgewidth=0.0,label='6T echo')
axarr[1,1].plot(intRange/sig,cumulIN/intRange,'b.-',label='file')
axarr[1,1].legend(loc='best')

plt.show()
