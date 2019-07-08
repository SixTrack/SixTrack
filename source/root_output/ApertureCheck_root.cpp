#include "ApertureCheck_root.h"

#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

ApertureCheckRootOutput* RootApertureOutput;

extern "C" void ApertureCheckRootInit()
{
    RootApertureOutput = new ApertureCheckRootOutput();
}

extern "C" void ApertureCheckWriteLossParticle(int turn_in, int i_in, int ix_in, char* bez_in, int bez_len, double slos_in, int ipart_in, double x_in, double xp_in, double y_in, double yp_in, double p_in, double dp_in, double ct_in, int a_in, int z_in, int q_in, int pdgid_in)
{
    RootApertureOutput->WriteLossParticle(turn_in, i_in, ix_in, bez_in, bez_len, slos_in, ipart_in, x_in, xp_in, y_in, yp_in, p_in, dp_in, ct_in, a_in, z_in, q_in, pdgid_in);
}

extern "C" void ApertureCheckWriteLossParticleF(int turn_in, int i_in, int ix_in, char* bez_in, int bez_len, double slos_in, int32_t fluka_uid_in, int32_t fluka_gen_in, double fluka_weight_in, double x_in, double xp_in, double y_in, double yp_in, double p_in, double dp_in, double ct_in, int a_in, int z_in, int q_in, int pdgid_in)
{
    RootApertureOutput->WriteLossParticleF(turn_in, i_in, ix_in, bez_in, bez_len, slos_in, fluka_uid_in, fluka_gen_in, fluka_weight_in, x_in, xp_in, y_in, yp_in, p_in, dp_in, ct_in, a_in, z_in, q_in, pdgid_in);
}

//Class functions

ApertureCheckRootOutput::ApertureCheckRootOutput()
{
    //Tree stuff
    ApertureLossTree = new TTree("ApertureLoss","ApertureLossTree");
    ApertureLossTree->Branch("turn",&turn,"turn/I");
    ApertureLossTree->Branch("i",&i,"i/I");
    ApertureLossTree->Branch("ix",&ix,"ix/I");
    ApertureLossTree->Branch("name",name,"name[49]/C");
    ApertureLossTree->Branch("ipart",&ipart,"ipart/I");
    ApertureLossTree->Branch("slos",&slos,"slos/D");
    ApertureLossTree->Branch("x", &x, "x/D");
    ApertureLossTree->Branch("xp",&xp,"xp/D");
    ApertureLossTree->Branch("y", &y, "y/D");
    ApertureLossTree->Branch("yp",&yp,"yp/D");

// This is the particle MOMENTUM in GeV
    ApertureLossTree->Branch("p",&p,"p/D");

// This is dp in eV
    ApertureLossTree->Branch("dp",&dp,"dp/D");

    ApertureLossTree->Branch("ct",&ct,"ct/D");
/*
S: a 16 bit signed integer
s: a 16 bit unsigned integer
I: a 32 bit signed integer
i: a 32 bit unsigned integer
L: a 64 bit signed integer
l: a 64 bit unsigned integer
*/
    ApertureLossTree->Branch("a",&na,"a/I");
    ApertureLossTree->Branch("z",&nz,"z/I");
    ApertureLossTree->Branch("q",&nq,"q/I");
    ApertureLossTree->Branch("pdgid",&pdgid,"pdgid/I");

//FLUKA variables
    ApertureLossTree->Branch("fluka_uid",&fluka_uid,"fluka_uid/i");
    ApertureLossTree->Branch("fluka_gen",&fluka_gen,"fluka_gen/i");
    ApertureLossTree->Branch("fluka_weight",&fluka_weight,"fluka_weight/D");
}

void ApertureCheckRootOutput::WriteLossParticle(int turn_in, int i_in, int ix_in, char* bez_in, int bez_len, double slos_in, int ipart_in, double x_in, double xp_in, double y_in, double yp_in, double p_in, double dp_in, double ct_in, int a_in, int z_in, int q_in, int pdgid_in)
{
    turn = turn_in;
    i = i_in;
    ix = ix_in;

    std::string tmpname(bez_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(name,tmpname.substr(0,bez_len).c_str(),49);

    slos = slos_in;
    ipart = ipart_in;
    x = x_in;
    xp = xp_in;
    y = y_in;
    yp = yp_in;
    p = p_in;
    dp = dp_in;
    ct = ct_in;

    na = a_in;
    nz = z_in;

    nq = q_in;
    pdgid = pdgid_in;

//    std::cout << "Writing Aperture loss: " << turn << " - " << slos << std::endl;
    //Do the write
    ApertureLossTree->Fill();
}

void ApertureCheckRootOutput::WriteLossParticleF(int turn_in, int i_in, int ix_in, char* bez_in, int bez_len, double slos_in, int32_t fluka_uid_in, int32_t fluka_gen_in, double fluka_weight_in, double x_in, double xp_in, double y_in, double yp_in, double p_in, double dp_in, double ct_in, int a_in, int z_in, int q_in, int pdgid_in)
{
    turn = turn_in;
    i = i_in;
    ix = ix_in;

    std::string tmpname(bez_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(name,tmpname.substr(0,bez_len).c_str(),49);

    slos = slos_in;

    fluka_uid = fluka_uid_in;
    fluka_gen = fluka_gen_in;
    fluka_weight = fluka_weight_in;

    x = x_in;
    xp = xp_in;
    y = y_in;
    yp = yp_in;
    p = p_in;
    dp = dp_in;
    ct = ct_in;

    na = a_in;
    nz = z_in;

    nq = q_in;
    pdgid = pdgid_in;

//    std::cout << "Writing Aperture loss: " << turn << " - " << slos << std::endl;
    //Do the write
    ApertureLossTree->Fill();
}


