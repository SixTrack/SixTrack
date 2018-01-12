#include "ApertureCheck_root.h"

#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

//Aperture check functions
Int_t turn;
Int_t i;
Int_t ix;
char loss_name[17];
Double_t slos;
Int_t ipart;
Double_t x;
Double_t xp;
Double_t y;
Double_t yp;
Double_t ct;
Double_t e;
Double_t dp;

TTree *ApertureLossTree;

extern "C" void ApertureCheckRootInit()
{
    //Tree stuff
    ApertureLossTree = new TTree("ApertureLoss","ApertureLossTree");
    ApertureLossTree->Branch("turn",&turn,"turn/I");
    ApertureLossTree->Branch("i",&i,"i/I");
    ApertureLossTree->Branch("ix",&ix,"ix/I");
    ApertureLossTree->Branch("name",&loss_name,"name[17]/C");
    ApertureLossTree->Branch("ipart",&ipart,"ipart/I");
    ApertureLossTree->Branch("slos",&slos,"slos/D");
    ApertureLossTree->Branch("x", &x, "x/D");
    ApertureLossTree->Branch("xp",&xp,"xp/D");
    ApertureLossTree->Branch("y", &y, "y/D");
    ApertureLossTree->Branch("yp",&yp,"yp/D");
    ApertureLossTree->Branch("e",&e,"e/D");
    ApertureLossTree->Branch("dp",&dp,"dp/D");
    ApertureLossTree->Branch("ct",&ct,"ct/D");
}

extern "C" void ApertureCheckWriteLossParticle(int turn_in, int i_in, int ix_in, char* bez_in, int bez_len, double slos_in, int ipart_in, double x_in, double xp_in, double y_in, double yp_in, double e_in, double dp_in, double ct_in)
{
    turn = turn_in;
    i = i_in;
    ix = ix_in;

    std::string tmpname(bez_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(loss_name,tmpname.substr(0,bez_len).c_str(),17);

    slos = slos_in;
    ipart = ipart_in;
    x = x_in;
    xp = xp_in;
    y = y_in;
    yp = yp_in;
    e = e_in;
    dp = dp_in;
    ct = ct_in;

//    std::cout << "Writing Aperture loss: " << turn << " - " << slos << std::endl;
    //Do the write
    ApertureLossTree->Fill();
}

