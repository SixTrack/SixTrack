import matplotlib.pyplot as plt
import numpy as np
import math

def fitProfile( x, fitParams, cLen ):
    '''
    fitParams[0]: const term
    fitParams[1]: linear term
    ...
    fitParams[-2]: n-th order term
    fitParams[-1]: scaling factor
    '''
    y=0.0
    for ii in range(len(fitParams)-1):
        if (ii==2):
            y+=(x**ii)*fitParams[ii]/cLen
        else:
            y+=(x**ii)*fitParams[ii]
    return y*fitParams[-1]

def rotateBy(xIn,yIn,skewAngle=0.0,direct=True):
    if (type(xIn) is list and type(yIn) is list):
        xOut=[]; yOut=[]
        for tmpX,tmpY in zip(xIn,yIn):
            xOut.append(tmpX*math.cos(skewAngle)+math.sin(skewAngle)*tmpY)
            yOut.append(tmpY*math.cos(skewAngle)-math.sin(skewAngle)*tmpX)
    else:
        xOut=xIn*math.cos(skewAngle)+math.sin(skewAngle)*yIn
        yOut=yIn*math.cos(skewAngle)-math.sin(skewAngle)*xIn
    return xOut,yOut
    
def parseFirstImpacts(iFileName='FirstImpacts.dat'):
    print 'parsing file %s ...'%(iFileName)
    data=[]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            if (line.startswith('#')): continue
            data.append([])
            tmpData=line.strip().split()
            for ii,tmpDatum in zip(range(len(tmpData)),tmpData):
                data[-1].append(float(tmpDatum))
                if ( ii<=3 ):
                    data[-1][-1]=int(data[-1][-1])
    print '...done - read %i lines.'%(len(data))
    return data

def parseFlukaImpacts(iFileName='FLUKA_impacts.dat'):
    print 'parsing file %s ...'%(iFileName)
    data=[]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            if (line.startswith('#')): continue
            data.append([])
            tmpData=line.strip().split()
            for ii,tmpDatum in zip(range(len(tmpData)),tmpData):
                data[-1].append(float(tmpDatum))
                if ( ii==0 or ii>=7 ):
                    data[-1][-1]=int(data[-1][-1])
    print '...done - read %i lines.'%(len(data))
    return data

def parseJawProfiles(iFileName='jaw_profiles.dat'):
    print 'parsing file %s ...'%(iFileName)
    data=[]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            if (line.startswith('#')): continue
            data.append([])
            tmpData=line.strip().split()
            for ii,tmpDatum in zip(range(len(tmpData)),tmpData):
                data[-1].append(float(tmpDatum))
                if ( ii<=2 ):
                    data[-1][-1]=int(data[-1][-1])
    print '...done - read %i lines.'%(len(data))
    return data

def getFittingData(iFileName='fort.6'):
    print 'parsing file %s ...'%(iFileName)
    profiles=[ [], [], [] ]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            line=line.strip()
            if (line.startswith('Fit point #')):
                data=line.split()
                profiles[0].append(float(data[4]))
                profiles[1].append(float(data[5])*1000)
                profiles[2].append(float(data[6])*1000)
    return profiles

def parseCollGaps(iFileName='collgaps.dat'):
    print 'parsing file %s ...'%(iFileName)
    collData=[]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            line=line.strip()
            if (line.startswith('#')): continue
            data=line.split()
            for ii in range(len(data)):
                if ( ii!=1 and ii !=6):
                    data[ii]=float(data[ii])
                    if ( ii==0 ):
                        data[ii]=int(data[ii]+1E-4)
            collData.append(data[:])
    print ' ...acquired %i collimators'%(len(collData))
    return collData

def parseCollDBtemp(iFileName='collimator-temp.db'):
    print 'parsing file %s ...'%(iFileName)
    collData=[]
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            line=line.strip()
            if (line.startswith('#')): 
                # a new collimator
                collData.append([])
            else:
                collData[-1].append(line)
                if (len(collData[-1])>2):
                    collData[-1][-1]=float(collData[-1][-1])
    print ' ...acquired %i collimators'%(len(collData))
    return collData

def plotFirstImpacts(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data):
    plt.figure('First impacts',figsize=(16,16))
    iPlot=0
    for iCol,lCol,jCol in zip(iCols,lCols,range(len(lCols))):
        Ss_in=[] ; Ss_out=[]
        Hs_in=[] ; Hs_out=[]
        Vs_in=[] ; Vs_out=[]
        for datum in data:
            if (datum[2]==iCol):
                # - s_in: entrance at jaw/slice or actual impact point;
                # - x_in,xp_in,y_in,yp_in (jaw ref sys): particle at front face of
                #   collimator/slice, not at impact on jaw
                # - x_out,xp_out,y_out,yp_out (jaw ref sys): particle exit point at
                #   collimator/slice, or coordinate of hard interaction
                ss=datum[4]
                xx=datum[6]
                yy=datum[8]
                Ss_in.append(ss)
                Hs_in.append(xx*1000)
                Vs_in.append(yy*1000)
                ss=datum[5]
                xx=datum[10]
                yy=datum[12]
                Ss_out.append(ss)
                Hs_out.append(xx*1000)
                Vs_out.append(yy*1000)
        
        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        if ( jCol==0 ):
            plt.plot(sixCurves[0],sixCurves[1],'ro-',label='6T fit-L')
            plt.plot(sixCurves[0],sixCurves[2],'mo-',label='6T fit-R')
        plt.plot(fitProfileS[jCol],fitProfileY1[jCol],'b-',label='expected-L')
        plt.plot(fitProfileS[jCol],fitProfileY2[jCol],'c-',label='expected-R')
        plt.plot(Ss_in,Hs_in,'go',label='in')
        plt.plot(Ss_out,Hs_out,'yo',label='out')
        plt.xlabel(r's_{jaw} [m]')
        plt.ylabel(r'x_{jaw} [mm]')
        plt.title('iColl=%i - Cleaning plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        plt.plot(Ss_in,Vs_in,'go',label='in')
        plt.plot(Ss_out,Vs_out,'yo',label='out')
        plt.xlabel(r's_{jaw} [m]')
        plt.ylabel(r'y_{jaw} [mm]')
        plt.title('iColl=%i - Ortoghonal plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        plt.plot(Hs_in,Vs_in,'go',label='in')
        plt.plot(Hs_out,Vs_out,'yo',label='out')
        plt.xlabel(r'x_{jaw} [mm]')
        plt.ylabel(r'y_{jaw} [mm]')
        plt.title('iColl=%i - transverse plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

    plt.show()

def plotFlukaImpacts(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data,skewAngles):
    plt.figure('FLUKA impacts',figsize=(16,16))
    iPlot=0
    for iCol,lCol,jCol,skewAngle in zip(iCols,lCols,range(len(lCols)),skewAngles):
        Ss=[]
        Hs=[]
        Vs=[]
        for datum in data:
            if (datum[0]==iCol):
                ss=datum[2]
                xx=datum[3]
                yy=datum[5]
                Ss.append(ss)
                Hs.append(xx)
                Vs.append(yy)
        
        if ( jCol==0 ):
            xx1,yy1=rotateBy(sixCurves[1],[0.0 for col in range(len(sixCurves[0]))],skewAngle=-skewAngle)
            xx2,yy2=rotateBy(sixCurves[2],[0.0 for col in range(len(sixCurves[0]))],skewAngle=-skewAngle)
        xp1,yp1=rotateBy(fitProfileY1[jCol],[0.0 for col in range(len(fitProfileS[jCol]))],skewAngle=-skewAngle)
        xp2,yp2=rotateBy(fitProfileY2[jCol],[0.0 for col in range(len(fitProfileS[jCol]))],skewAngle=-skewAngle)
            
        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        if ( jCol==0 ):
            plt.plot(sixCurves[0],xx1,'ro-',label='6T fit-L')
            plt.plot(sixCurves[0],xx2,'mo-',label='6T fit-R')
        plt.plot(Ss,Hs,'go',label='impact')
        plt.plot(fitProfileS[jCol],xp1,'b-',label='expected-L')
        plt.plot(fitProfileS[jCol],xp2,'c-',label='expected-R')
        plt.xlabel(r's [m]')
        plt.ylabel(r'x [mm]')
        plt.title('iColl=%i - x-s plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        if ( jCol==0 ):
            plt.plot(sixCurves[0],yy1,'ro-',label='6T fit-L')
            plt.plot(sixCurves[0],yy2,'mo-',label='6T fit-R')
        plt.plot(Ss,Vs,'go',label='impact')
        plt.plot(fitProfileS[jCol],yp1,'b-',label='expected-L')
        plt.plot(fitProfileS[jCol],yp2,'c-',label='expected-R')
        plt.xlabel(r's [m]')
        plt.ylabel(r'y [mm]')
        plt.title('iColl=%i - y-s plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        if ( jCol==0 ):
            plt.plot(xx1,yy1,'ro-',label='6T fit-L')
            plt.plot(xx2,yy2,'mo-',label='6T fit-R')
        plt.plot(Hs,Vs,'go',label='impact')
        plt.plot(xp1,yp1,'b-',label='expected-L')
        plt.plot(xp2,yp2,'c-',label='expected-R')
        plt.xlabel(r'x [mm]')
        plt.ylabel(r'y [mm]')
        plt.title('iColl=%i - x-y plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

    plt.show()

def plotJawProfiles(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data):
    plt.figure('Jaw Profiles',figsize=(16,16))
    iPlot=0
    for iCol,lCol,jCol in zip(iCols,lCols,range(len(lCols))):
        Ss_in=[] ; Ss_out=[]
        Hs_in=[] ; Hs_out=[]
        Vs_in=[] ; Vs_out=[]
        for datum in data:
            if (datum[0]==iCol):
                ss=datum[7]
                xx=datum[3]
                yy=datum[5]
                if ( datum[-1]==1 ):
                    # entrance
                    Ss_in.append(ss)
                    Hs_in.append(xx*1000)
                    Vs_in.append(yy*1000)
                else:
                    # exit
                    Ss_out.append(ss)
                    Hs_out.append(xx*1000)
                    Vs_out.append(yy*1000)
        
        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        if ( jCol==0 ):
            plt.plot(sixCurves[0],sixCurves[1],'ro-',label='6T fit-L')
            plt.plot(sixCurves[0],sixCurves[2],'mo-',label='6T fit-R')
        plt.plot(fitProfileS[jCol],fitProfileY1[jCol],'b-',label='expected-L')
        plt.plot(fitProfileS[jCol],fitProfileY2[jCol],'c-',label='expected-R')
        plt.plot(Ss_in,Hs_in,'go',label='in')
        plt.plot(Ss_out,Hs_out,'yo',label='out')
        plt.xlabel(r's_{jaw} [m]')
        plt.ylabel(r'x_{jaw} [mm]')
        plt.title('iColl=%i - Cleaning plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        plt.plot(Ss_in,Vs_in,'go',label='in')
        plt.plot(Ss_out,Vs_out,'yo',label='out')
        plt.xlabel(r's_{jaw} [m]')
        plt.ylabel(r'y_{jaw} [mm]')
        plt.title('iColl=%i - Ortoghonal plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

        iPlot+=1
        plt.subplot(len(iCols),3,iPlot)
        plt.plot(Hs_in,Vs_in,'go',label='in')
        plt.plot(Hs_out,Vs_out,'yo',label='out')
        plt.xlabel(r'x_{jaw} [mm]')
        plt.ylabel(r'y_{jaw} [mm]')
        plt.title('iColl=%i - transverse plane'%(iCol))
        plt.legend(loc='best',fontsize=10)
        plt.tight_layout()
        plt.grid()

    plt.show()

if ( __name__ == "__main__" ):
    iCols=[9,10,11,15]
    fitParams=[
        # iCol=9
        [   2.70E-3, -0.18,  0.18,  0.0, 0.0, 0.0, 2 ], 
        [ 5.1962E-4,  0.09, -0.27, 50.0, 0.0, 0.0, 2 ],
        # iCol=10
        [ 0., 1. ],
        [ 0., 1. ],
        # iCol=11
        [ 0., 1. ],
        [ 0., 1. ],
        # iCol=15
        [ 0., 1. ],
        [ 0., 1. ]
    ]
    nPoints=110
    lDebug=False
    
    collData=parseCollGaps()
    lCols=[] ; hGaps=[] ; skewAngles=[] ; tilt1=[] ; tilt2=[]
    for iCol in iCols:
        for collDatum in collData:
            if ( collDatum[0]==iCol ):
                lCols.append(collDatum[7])
                hGaps.append(collDatum[5])
                skewAngles.append(collDatum[2])
                tilt1.append(collDatum[10])
                tilt2.append(collDatum[11])
                break

    collDBdata=parseCollDBtemp()
    offSets=[]
    for iCol in iCols:
        offSets.append(collDBdata[iCol-1][4])

    if (lDebug):
        for iCol,jCol in zip(iCols,range(len(iCols))):
            print iCol, lCols[jCol], hGaps[jCol], skewAngles[jCol], tilt1[jCol], tilt2[jCol], offSets[jCol]
            
    fitProfileS=[] ; fitProfileY1=[] ; fitProfileY2=[]
    for lCol,hGap,jCol in zip(lCols,hGaps,range(len(iCols))):
        Ss=[ float(ii)/nPoints*lCol for ii in range(nPoints+1) ]
        if ( tilt1[jCol]>0 ):
            Y1s=[ ( fitProfile(s,fitParams[  2*jCol],lCol)+hGap+tilt1[jCol]*(s     )+offSets[jCol])*1000 for s in Ss ]
        else:
            Y1s=[ ( fitProfile(s,fitParams[  2*jCol],lCol)+hGap+tilt1[jCol]*(s-lCol)+offSets[jCol])*1000 for s in Ss ]
        if ( tilt2[jCol]>0 ):
            Y2s=[ (-fitProfile(s,fitParams[1+2*jCol],lCol)-hGap+tilt2[jCol]*(s-lCol)+offSets[jCol])*1000 for s in Ss ]
        else:
            Y2s=[ (-fitProfile(s,fitParams[1+2*jCol],lCol)-hGap+tilt2[jCol]*(s     )+offSets[jCol])*1000 for s in Ss ]
        fitProfileS.append(Ss[:])
        fitProfileY1.append(Y1s[:])
        fitProfileY2.append(Y2s[:])

    sixCurves=getFittingData()
    sixCurves[1]=[tmpPoint+offSets[0]*1000 for tmpPoint in sixCurves[1]]
    sixCurves[2]=[tmpPoint+offSets[0]*1000 for tmpPoint in sixCurves[2]]
    
    data=parseFirstImpacts()
    plotFirstImpacts(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data)
    data=parseFlukaImpacts()
    plotFlukaImpacts(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data,skewAngles)
    data=parseJawProfiles()
    plotJawProfiles(iCols,lCols,fitProfileS,fitProfileY1,fitProfileY2,sixCurves,data)
