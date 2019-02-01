#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#plotType = None
#plotType="IPACsingleCol_crabVoltage"
#plotType="IPACsingleCol_crabFail2"
#plotType="manual_crabfail"
plotType="generic"

if plotType=="IPACsingleCol_crabVoltage":
    textwidth = 3.25 #inches, for 2Dpic paper
    
    from matplotlib import rcParams,rc
#    rcParams.update({'text.usetex': True}) #slow
    DPI = 300
    rc('font',**{'family':'serif','serif':['Times'],'size':10})
    rcParams['figure.figsize'] = textwidth, textwidth/1.618

elif plotType=="IPACsingleCol_crabFail1":
    textwidth = 3.25 #inches, for 2Dpic paper
    
    from matplotlib import rcParams,rc
#    rcParams.update({'text.usetex': True}) #slow
    DPI = 300
    rc('font',**{'family':'serif','serif':['Times'],'size':10})
    rcParams['figure.figsize'] = textwidth, textwidth/1.618

    (fig_fail, axes_fail) = plt.subplots(nrows=2, ncols=1,figsize=(textwidth,textwidth/1.618*2), sharex=True, sharey=False, num=1000,dpi=DPI)

elif plotType=="IPACsingleCol_crabFail2":
#    textwidth = 3.25 #inches, for 2Dpic paper
    textwidth = 483.69684/72.27
    
    from matplotlib import rcParams,rc
#    rcParams.update({'text.usetex': True}) #slow
    DPI = 100
    rc('font',**{'family':'serif','serif':['Times'],'size':10})
    rcParams['figure.figsize'] = textwidth, textwidth/1.618
    
    turn_volt = None
    turn_phase = None
    volt_volt = None
    phase_phase = None

elif plotType=="manual_crabfail":
    textwidth = 483.69684/72.27
    from matplotlib import rcParams,rc
#    rcParams.update({'text.usetex': True}) #slow
    DPI = 300
    rc('font',**{'family':'serif','serif':['Times'],'size':10})
    rcParams['figure.figsize'] = textwidth, textwidth/1.618

else:
    DPI=300

class DYNKdata:
    #Data for one (element,attribute)
    
    element   = None
    attribute = None
    
    turn    = None # Array of turn numbers
    setIDX  = None # Array of setIDX
    funName = None # Array of funname
    data    = None # Matrix of data
    #How many elements are we setting = how many columns in data
    datalen = None 

    def __init__(self, element, attribute,datalen):
        self.element   = element
        self.attribute = attribute

        self.turn    = []
        self.setIDX  = []
        self.funName = []
        self.data    = [] #First index: turn, second, element instance
        
        self.datalen = datalen
        
    def toNumarray(self):
        self.turn = np.asarray(self.turn)
        self.data = np.asarray(self.data)
        assert self.datalen == self.data.shape[1]
    
dynkData = {}

ifile = open("dynksets.dat", 'r')
lnum = 1
for l in ifile:
    l = l.strip()
    if l[0] == "#":
        continue
    ls = l.split()
    
    TURN = int(ls[0])
    EL_ATT = (ls[1],ls[2])
    SETIDX = int(ls[3])
    FUNCTION = ls[4]
    if len(ls) == 6:
        if lnum ==1:
            print "New file format"
        DATALEN = 1
        DATA = [float(ls[5])]
    elif len(ls) >= 7: 
        if lnum ==1:
            print "Old file format"
        DATALEN = int(ls[5])
        DATA = map(float,ls[6:])
    else:
        print "ERROR: Invalid file format"
        exit(1)
    
    if DATALEN != len(DATA):
        print "ERROR! DATALEN != DATA"
        exit(1)

    if not EL_ATT in dynkData:
        dynkData[EL_ATT] = DYNKdata(EL_ATT[0],EL_ATT[1],DATALEN)

    dynkData[EL_ATT].turn.append(TURN)
    dynkData[EL_ATT].setIDX.append(SETIDX)
    dynkData[EL_ATT].funName.append(FUNCTION)
    dynkData[EL_ATT].data.append(DATA)

    lnum = lnum+1

ifile.close()

if len(dynkData) == 0:
    print "No data"
    exit(0)

for (k,i) in zip(dynkData,xrange(len(dynkData))):
    dynkData[k].toNumarray()
    
    #Separate figure:
    plt.figure(10+i)
    plt.plot(dynkData[k].turn,dynkData[k].data[:,0])
    plt.xlabel("Turn")
    plt.ylabel("Setting")
    plt.title("%s:%s" %(k[0],k[1]))
    #  plt.grid(True)
    # plt.xticks(np.arange(dynkData[k].turn[0], dynkData[k].turn[-1], 1.0))
    plt.savefig("dynksets-%s_%s.png" %(k[0],k[1]),dpi=DPI)
    
    if plotType=="IPACsingleCol_crabFail1":
        if k[0] == "CRAB_IP1_L1" and k[1] == "voltage":
            axes_fail[0].plot(dynkData[k].turn,dynkData[k].data[:,0])
            axes_fail[0].set_ylabel("Voltage [MV]")
            axes_fail[0].yaxis.set_label_coords(-0.15,0.5)

        if k[0] == "CRAB_IP1_L1" and k[1] == "phase":
            axes_fail[1].plot(dynkData[k].turn,dynkData[k].data[:,0])
            axes_fail[1].set_xlabel("Turn")
            axes_fail[1].set_ylabel("Phase [radians]")
#            axes_fail[0].yticks([0,np.pi/2,2*np.pi/2,3*np.pi/2,4*np.pi/2,5*np.pi/2,6*np.pi/2],[r"$0$",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$",r"$5\pi/2$",r"$6\pi/2$"])
            axes_fail[1].yaxis.set_ticks([0,np.pi/2,2*np.pi/2,3*np.pi/2,4*np.pi/2,5*np.pi/2,6*np.pi/2])
            axes_fail[1].yaxis.set_ticklabels([r"$0$",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$",r"$5\pi/2$",r""])
            axes_fail[1].yaxis.set_label_coords(-0.15,0.5)
    elif plotType=="IPACsingleCol_crabFail2":
        if k[0] == "CRAB_IP1_L1" and k[1] == "voltage":
            turn_volt = dynkData[k].turn
            volt_volt = dynkData[k].data[:,0]

        elif k[0] == "CRAB_IP1_L1" and k[1] == "phase":
            turn_phase  = dynkData[k].turn
            phase_phase = dynkData[k].data[:,0]

    #Combined figure
    plt.figure(1)
    if plotType == "IPACsingleCol_crabVoltage":
        plotLabel = "%s"%(k[0])
        plt.plot(dynkData[k].turn,dynkData[k].data[:,0], label=plotLabel[-2:] )
    else:
        plotLabel = "%s:%s"%(k[0],k[1])
        print plotLabel
        plt.plot(dynkData[k].turn,dynkData[k].data[:,0], label=plotLabel )
    #plt.xticks(np.arange(dynkData[k].turn[0], dynkData[k].turn[-1], 1.0))

#Final touches
plt.figure(1)
plt.xlabel("Turn")
if plotType == "IPACsingleCol_crabVoltage":
    plt.ylabel("Voltage [MV]")

    handles, labels = plt.gca().get_legend_handles_labels()
    # sort them by labels
    import operator
    hl = sorted(zip(handles, labels),
                key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)

    plt.legend(handles2, labels2, loc=4,frameon=False,ncol=2,fontsize='small')

    plt.subplots_adjust(left=0.16,bottom=0.18, right=0.95,top=0.95)
    plt.savefig("dynksets_all.png",dpi=DPI)
else:
    plt.ylabel("Setting")
    #plt.legend(loc=0,frameon=False,ncol=2)
    plt.legend(loc=0,frameon=False)

    plt.savefig("dynksets_all.png",dpi=DPI)
plt.grid(True)

if plotType == "IPACsingleCol_crabFail1":
    fig_fail.subplots_adjust(hspace=0,left=0.17,right=0.96,top=0.98,bottom=0.1)
    fig_fail.savefig("fail_voltagePhase.png",dpi=DPI)

elif plotType=="IPACsingleCol_crabFail2":
    plt.figure(1000,dpi=DPI)
    ax1=plt.gca()
    ax2=plt.twinx()

    p1, = ax1.plot(turn_volt, volt_volt, 'b-',label="voltage")
    ax1.set_ylabel("Voltage [MV]")
    ax1.set_xlabel("Turn")

    p2, = ax2.plot(turn_phase, phase_phase,'r--',label="Phase")
    ax2.set_ylabel("Phase [radians]",color='r')
    ax2.yaxis.set_ticks([0,np.pi/2,2*np.pi/2,3*np.pi/2,4*np.pi/2,5*np.pi/2,6*np.pi/2])
    ax2.yaxis.set_ticklabels([r"$0$",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$",r"$5\pi/2$",r"$3\pi$"])

    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    ax2.legend((p1,p2),("Voltage","Phase"),loc=2,frameon=False,fontsize='small')

    plt.subplots_adjust(left=0.15,bottom=0.175, right=0.83,top=0.96)

    plt.savefig("fail_voltagePhase2.png",dpi=DPI)
    plt.savefig("fail_voltagePhase2.eps")

plt.show()
