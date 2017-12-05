#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

//Optics functions
Int_t n;
Double_t s;
Double_t orbit_x;
Double_t orbit_y;
Double_t orbit_xp;
Double_t orbit_yp;
Double_t beta_x;
Double_t beta_y;
Double_t alpha_x;
Double_t alpha_y;
//Double_t gamma_x;
//Double_t gamma_y;
Double_t dispersion_x;
Double_t dispersion_y;

Double_t dispersion_xp;
Double_t dispersion_yp;

//Double_t mu_x;
//Double_t mu_y;

Char_t name[17];

TTree *OpticsTree;

//Dumps optical functions
extern "C" void OpticsRootInit()
{
    //Tree stuff
    OpticsTree = new TTree("LinearOptics","LinearOpticsTree");
    OpticsTree->Branch("n",&n,"n/I");

    OpticsTree->Branch("s", &s, "s/D");

    OpticsTree->Branch("orbit_x", &orbit_x, "orbit_x/D");
    OpticsTree->Branch("orbit_y", &orbit_y, "orbit_y/D");

    OpticsTree->Branch("orbit_xp", &orbit_xp, "orbit_xp/D");
    OpticsTree->Branch("orbit_yp", &orbit_yp, "orbit_yp/D");

    OpticsTree->Branch("beta_x", &beta_x, "beta_x/D");
    OpticsTree->Branch("beta_y", &beta_y, "beta_y/D");

    OpticsTree->Branch("alpha_x", &alpha_x, "alpha_x/D");
    OpticsTree->Branch("alpha_y", &alpha_y, "alpha_y/D");

    OpticsTree->Branch("dispersion_x", &dispersion_x, "dispersion_x/D");
    OpticsTree->Branch("dispersion_y", &dispersion_y, "dispersion_y/D");

    OpticsTree->Branch("dispersion_xp", &dispersion_xp, "dispersion_xp/D");
    OpticsTree->Branch("dispersion_yp", &dispersion_yp, "dispersion_yp/D");

    OpticsTree->Branch("name", name, "name[17]/C");
//    OpticsTree->Branch("mu_x", &mu_x, "mu_x/D");
//    OpticsTree->Branch("mu_y", &mu_y, "mu_y/D");
}

extern "C" void OpticsRootWrite(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in)
{
    n = n_in;
    std::string tmpname(name_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(name,tmpname.substr(0,c_len).c_str(),17);

    s = s_in;

    orbit_x  = x_in;
    orbit_y  = y_in;

    orbit_xp = xp_in;
    orbit_yp = yp_in;

    beta_x = beta_x_in;
    beta_y = beta_y_in;

    alpha_x = alpha_x_in;
    alpha_y = alpha_y_in;

    dispersion_x = dispersion_x_in;
    dispersion_y = dispersion_y_in;

    dispersion_xp = dispersion_xp_in;
    dispersion_yp = dispersion_yp_in;

    //Do the write
    OpticsTree->Fill();
}

