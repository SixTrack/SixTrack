
class test():
    def __init__(self):
        self.name=''
        self.duration=0 # [s]
    @staticmethod
    def fromLineCTestCostData(line): # acquire a single line
        data=line.strip().split()
        newTest=test()
        newTest.name=data[0]
        newTest.duration=float(data[2])
        return newTest
    @staticmethod
    def fromLineCTestSTDOUT(line): # acquire a single line
        data=line.strip().split()
        newTest=test()
        newTest.name=data[3]
        newTest.duration=float(data[6])
        return newTest

def runTests(cmd='ctest -L fast -j 6',nTimes=5,dest='../..'):
    '''
    run the tests and copy the CTestCostData.txt/STDOUT in a specific folder
    ctest is then run nTimes, and each file will be saved with a time stamp, eg:
    commit_5d495f366e641e0c7ac972b9a422f74d296d02e6/
    |_ ctest_STDOUT_2018-11-14_15-40-32.txt
    |_ ctest_STDOUT_2018-11-14_15-42-50.txt
    ...

    NB: it is not at all clear to me why:
    subprocess.call(['ctest','-L fast','-j 6'])
    leads to executing all tests and not only the fast ones...
    '''
    import subprocess
    import datetime
    import os
    if (nTimes>0):
        if ( not os.path.isdir(dest) ):
            os.makedirs(dest)
        for ii in range(nTimes):
            now=datetime.datetime.now()
            print ' running %s at %s - iteration # %i'%(cmd,now.strftime('%Y-%m-%d %H:%M:%S'),ii+1)
            oFileName='%s/ctest_STDOUT_%s.txt'%(dest,now.strftime('%Y-%m-%d_%H-%M-%S'))
            subprocess.call('%s > %s'%(cmd,oFileName),shell=True)
            #subprocess.call('cp Testing/Temporary/CTestCostData.txt %s/CTestCostData_%s.txt'%(dest,now.strftime('%Y-%m-%d_%H-%M-%S')),shell=True)

def acquireDataSetCostData(files='../CTestCostData_*.txt'):
    '''
    parse CTestCostData.txt files
    returns testSet[testName][fileName]
    '''
    import glob
    testSet={}
    for fname in glob.glob(files): # parse all available files
        tmpTestSet=parseCTestCostDataFile(fname)
        for tmpTest in tmpTestSet.values():
            if ( not testSet.has_key(tmpTest.name) ):
                testSet[tmpTest.name]={}
            testSet[tmpTest.name][fname]=tmpTest
    return testSet
        
def acquireDataSetCTestSTDOUT(files='../ctest_STDOUT_*.txt'):
    '''
    parse STDOUT files of ctest
    returns testSet[testName][fileName]
    '''
    import glob
    testSet={}
    for fname in glob.glob(files): # parse all available files
        tmpTestSet=parseCTestSTDOUTFile(fname)
        for tmpTest in tmpTestSet.values():
            if ( not testSet.has_key(tmpTest.name) ):
                testSet[tmpTest.name]={}
            testSet[tmpTest.name][fname]=tmpTest
    return testSet
        
def parseCTestCostDataFile(iFileName):
    '''
    parse a single CTestCostData.txt file
    returns testSet[testName]
    '''
    print 'getting ctest times from file %s ...'%(iFileName)
    tests={}
    nFailed=0
    with open(iFileName,'r') as iFile:
        lFailed=False
        for line in iFile.readlines():
            if (line.startswith('---')):
                lFailed=True
                continue
            if ( lFailed ):
                nFailed+=1
            else:
                tmpTest=test.fromLineCTestCostData(line)
                tests[tmpTest.name]=tmpTest
    print 'done - acquired %i tests, out of which %i failed.'%(len(tests),nFailed)
    return tests

def parseCTestSTDOUTFile(iFileName):
    '''
    parse a single STDOUT file of ctest
    returns testSet[testName]
    '''
    print 'getting ctest times from file %s ...'%(iFileName)
    tests={}
    nFailed=0
    with open(iFileName,'r') as iFile:
        for line in iFile.readlines():
            if ('Passed' in line):
                tmpTest=test.fromLineCTestSTDOUT(line)
                tests[tmpTest.name]=tmpTest
            elif ('Failed' in line):
                nFailed+=1
    print 'done - acquired %i tests, out of which %i failed.'%(len(tests),nFailed)
    return tests

def buildStats(testSet,testNames=[],lSum=False):
    import math
    stats={}
    for testName in testNames:
        print 'computing statistics for test %s ...'%(testName)
        stats[testName]={}
        durations=[ testSet[testName][fname].duration for fname in testSet[testName].keys() ]
        stats[testName]['N']=len(durations)
        stats[testName]['min']=min(durations)
        stats[testName]['max']=max(durations)
        stats[testName]['ave']=sum(durations)/stats[testName]['N']
        stats[testName]['err']=math.sqrt(sum([duration**2 for duration in durations])/stats[testName]['N']-stats[testName]['ave']**2)
    if lSum:
        print 'computing statistics on sum of all tests...'
        stats['sum']={}
        stats['sum']['min']=sum([ stats[testName]['min'] for testName in testNames])
        stats['sum']['max']=sum([ stats[testName]['max'] for testName in testNames])
        stats['sum']['ave']=sum([ stats[testName]['ave'] for testName in testNames])
        stats['sum']['err']=math.sqrt(sum([ stats[testName]['err']**2 for testName in testNames]))
        stats['sum']['N']=len(testNames)
    return stats

def plotSingleTests(stats,testNames=[],commitHashes=[],commitNames=[],oFileName='output.png'):
    import matplotlib.pyplot as plt
    f, axarr = plt.subplots(1,1)
    plt.rcParams.update({'font.size': 10})
    # axarr.set_title(colliName)
    axarr.set_xlabel('commits')
    axarr.set_ylabel('duration [s]')
    axarr.grid()
    NUM_COLORS=len(testNames)
    cm = plt.get_cmap('gist_rainbow')
    axarr.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    Xs=[ii for ii in range(len(commitNames))]
    axarr.set_xlim(Xs[0]-0.5,Xs[-1]+0.5)
    plt.xticks(Xs,commitNames,rotation=90)
    for testName in testNames:
        y=[stats[commitHash][testName]['ave'] for commitHash in commitHashes if stats[commitHash].has_key(testName) ]
        yerr=[stats[commitHash][testName]['err'] for commitHash in commitHashes if stats[commitHash].has_key(testName) ]
        x=[ tmpX for tmpX,commitHash in zip(Xs,commitHashes) if stats[commitHash].has_key(testName) ]
        if (len(x)>0 ):
            axarr.errorbar(x,y,yerr=yerr,label='%s - ave - %i points'%(testName,stats[commitHash][testName]['N']),fmt='o')
    for testName in testNames:
        y=[stats[commitHash][testName]['min'] for commitHash in commitHashes if stats[commitHash].has_key(testName) ]
        x=[ tmpX for tmpX,commitHash in zip(Xs,commitHashes) if stats[commitHash].has_key(testName) ]
        if (len(x)>0 ):
            axarr.plot(x,y,marker='*')
    for testName in testNames:
        y=[stats[commitHash][testName]['max'] for commitHash in commitHashes if stats[commitHash].has_key(testName) ]
        x=[ tmpX for tmpX,commitHash in zip(Xs,commitHashes) if stats[commitHash].has_key(testName) ]
        if (len(x)>0 ):
            axarr.plot(x,y,marker='*')
    axarr.legend(loc='best')
    plt.savefig(oFileName,bbox_inches='tight')
    plt.show()
    plt.close()
        
if ( __name__ == '__main__' ):
    import collections
    # ctest command and number of repetitions:
    cmd='ctest -L fast -j 6'
    nTimes=5
    
    # what to plot
    # - commits (example):
    commmitHashes=collections.OrderedDict()
    commmitHashes['master, 2018-11-09']='5d495f366e641e0c7ac972b9a422f74d296d02e6'
    commmitHashes['refactoring code']='cecf082ef5d6c7d28be4f8c7176afa2b3ba428a7'
    commmitHashes['DA aper check in include']='cc9a2a233716e83cc092455314c1254c6325da00'
    commmitHashes['moving do_coll']='b612757c280b4e3bbce273fa557cf5ed7b0fe1b7'
    commmitHashes['separated llost']='e0c333fa2f798ced9c3d1aea24a6c620001c1e81'
    commmitHashes['interface aperture_reportLoss']='674e91b7a0f95eafadeac325d645d0107ac7974e'
    commmitHashes['restoring llostp']='7dca82d467199d9b0b97d8628aa851fca36a919b'
    lastCommitHash=commmitHashes['restoring llostp']
    # - which tests (can contain 'all')
    testNames=[
        'thin6d_ions',
        'thin4d_ions',
        'lost',        # with LIMI block
        'lostnumxv',    # with LIMI block
        'scatter_aperture'  # with LIMI block
        # 'orbit6d-ions-long', # 20M turns (non-fast)
        # 'prob1',  # 1M turns (non-fast)
        # 'prob3',  # 1M turns (non-fast)
    ]
    # plot also sum of all tests
    lSum=True
    
    # run tests
    if ( lastCommitHash is not None ):
        runTests(cmd=cmd,nTimes=nTimes,dest='../../commit_%s'%(lastCommitHash))
        if (lastCommitHash not in commmitHashes.values()):
            commmitHashes['latest commit']=lastCommitHash

    # compute averages and std deviations
    stats={}
    for commitLabel,commitHash in commmitHashes.iteritems():
        #testSet=acquireDataSetCostData(files='../../commit_%s/CTestCostData_*.txt'%(commitHash))
        testSet=acquireDataSetCTestSTDOUT(files='../../commit_%s/ctest_STDOUT_*.txt'%(commitHash))
        stats[commitHash]=buildStats(testSet,testNames=testSet.keys(),lSum=lSum)

    # plot
    if ( 'all' in testNames ):
        testNames=testSet.keys()
    # subset=[]
    # iSub=0
    # for testName in testNames:
    #     subset.append(testName)
    #     if (len(subset)==1):
    #         iSub+=1
    #         plotSingleTests(stats,testNames=subset,commitHashes=commmitHashes.values(),commitNames=commmitHashes.keys(),oFileName='set_%i.png'%(iSub))
    #         subset=[]
    # if (len(subset)!=0):
    #     iSub+=1
    #     plotSingleTests(stats,testNames=subset,commitHashes=commmitHashes.values(),commitNames=commmitHashes.keys(),oFileName='set_%i.png'%(iSub))
    #     subset=[]
    for testName in testNames:
        plotSingleTests(stats,testNames=[testName],commitHashes=commmitHashes.values(),commitNames=commmitHashes.keys(),oFileName='%s.png'%(testName))
    if (lSum):
        plotSingleTests(stats,testNames=['sum'],commitHashes=commmitHashes.values(),commitNames=commmitHashes.keys(),oFileName='sum.png')
