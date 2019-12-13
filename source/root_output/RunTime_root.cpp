#include "TTree.h"
#include "TROOT.h"
#include "RunTime_root.h"

double pretime;
double trtime;
double posttime;
double totaltime;

TTree *RunTimeTree;

//Dumps optical functions
extern "C" void RunTimeRootInit()
{
    //Tree stuff
    RunTimeTree = new TTree("RunTime","RunTimeTree");
    RunTimeTree->Branch("pretime",&pretime,"pretime/D");
    RunTimeTree->Branch("trtime",&trtime,"trtime/D");
    RunTimeTree->Branch("posttime",&posttime,"posttime/D");
    RunTimeTree->Branch("totaltime",&totaltime,"totaltime/D");
}

extern "C" void RunTimeRootWrite(double pretime_in, double trtime_in, double posttime_in)
{
    pretime = pretime_in;
    trtime = trtime_in;
    posttime = posttime_in;
    totaltime = pretime_in + trtime_in + posttime_in;

    //Do the write
    RunTimeTree->Fill();
}

