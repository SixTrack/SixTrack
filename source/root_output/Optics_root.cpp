#include "TTree.h"
#include "TROOT.h"

#include "Optics_root.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

OpticsRootOutput* RootOptics; 

//Dumps optical functions
extern "C" void OpticsRootInit()
{
    RootOptics = new OpticsRootOutput();
}

extern "C" void OpticsRootWriteLin(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in)
{
    RootOptics->WriteLin(n_in, name_in, c_len, s_in, x_in, xp_in, y_in, yp_in, beta_x_in, beta_y_in, alpha_x_in, alpha_y_in, dispersion_x_in, dispersion_y_in, dispersion_xp_in, dispersion_yp_in);
}

//phi,bexi,bexii,bezi,bezii, alxi,alxii,alzi, &
//     &alzii, gaxi,gaxii,gazi,gazii,phxi,phxii,phzi,phzii, phxpi,        &
//     &phxpii,phzpi,phzpii,couuang,t(6,1),t(6,2),t(6,3),t(6,4),t(1,1),   &
//     &t(1,2),t(1,3),t(1,4)
extern "C" void OpticsRootWriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alxi,double alxii,double alzi, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double phxpii, double phzpi, double phzpii,  double couuang_in, double t61_in, double t62_in, double t63_in, double t64_in, double t11_in, double t12_in, double t13_in, double t14_in)
{
    RootOptics->WriteCpl(phi1, phi2, bexi, bexii, bezi, bezii, alxi, alxii, alzi, alzii, gaxi, gaxii, gazi, gazii, phxi, phxii, phzi, phzii, phxpi, phxpii, phzpi, phzpii,  couuang_in, t61_in, t62_in, t63_in, t64_in, t11_in, t12_in, t13_in, t14_in);
}

extern "C" void OpticsRootWrite()
{
    RootOptics->Fill();
}

//Class functions
OpticsRootOutput::OpticsRootOutput()
{
    //Tree stuff
    OpticsTree = new TTree("LinearOptics","LinearOpticsTree");
    OpticsTree->Branch("n",&n,"n/I");

    OpticsTree->Branch("name", name, "name[49]/C");

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


    //These all should be init to 0
    OpticsTree->Branch("mu_x", &mu_x, "mu_x/D");
    OpticsTree->Branch("mu_y", &mu_y, "mu_y/D");

    OpticsTree->Branch("beta_c_x", &beta_c_x, "beta_c_x/D");
    OpticsTree->Branch("beta_c_y", &beta_c_y, "beta_c_y/D");
    OpticsTree->Branch("beta_c_z", &beta_c_z, "beta_c_z/D");
    OpticsTree->Branch("beta_c_dp", &beta_c_dp, "beta_c_dp/D");

    OpticsTree->Branch("alpha_c_x", &alpha_c_x, "alpha_c_x/D");
    OpticsTree->Branch("alpha_c_y", &alpha_c_y, "alpha_c_y/D");
    OpticsTree->Branch("alpha_c_z", &alpha_c_z, "alpha_c_z/D");
    OpticsTree->Branch("alpha_c_dp", &alpha_c_dp, "alpha_c_dp/D");

    OpticsTree->Branch("gamma_c_x", &gamma_c_x, "gamma_c_x/D");
    OpticsTree->Branch("gamma_c_y", &gamma_c_y, "gamma_c_y/D");
    OpticsTree->Branch("gamma_c_z", &gamma_c_z, "gamma_c_z/D");
    OpticsTree->Branch("gamma_c_dp", &gamma_c_dp, "gamma_c_dp/D");

    OpticsTree->Branch("phi_x", &phi_x, "phi_x/D");
    OpticsTree->Branch("phi_y", &phi_y, "phi_y/D");
    OpticsTree->Branch("phi_z", &phi_z, "phi_z/D");
    OpticsTree->Branch("phi_dp", &phi_dp, "phi_dp/D");

    OpticsTree->Branch("phi_xp", &phi_xp, "phi_xp/D");
    OpticsTree->Branch("phi_yp", &phi_yp, "phi_yp/D");
    OpticsTree->Branch("phi_zp", &phi_zp, "phi_zp/D");
    OpticsTree->Branch("phi_dpp", &phi_dpp, "phi_dpp/D");

    OpticsTree->Branch("couuang", &couuang, "couuang/D");

    OpticsTree->Branch("t11", &t11, "t11/D");
    OpticsTree->Branch("t12", &t12, "t12/D");
    OpticsTree->Branch("t13", &t13, "t13/D");
    OpticsTree->Branch("t14", &t14, "t14/D");

    OpticsTree->Branch("t61", &t61, "t61/D");
    OpticsTree->Branch("t62", &t62, "t62/D");
    OpticsTree->Branch("t63", &t63, "t63/D");
    OpticsTree->Branch("t64", &t64, "t64/D");

    phi_x   = 0;
    phi_y   = 0;
    phi_z   = 0;
    phi_dp = 0;

    phi_xp  = 0;
    phi_yp  = 0;
    phi_zp  = 0;
    phi_dpp = 0;

    mu_x = 0;
    mu_y = 0;

    beta_c_x  = 0;
    beta_c_y  = 0;
    beta_c_z  = 0;
    beta_c_dp = 0;

    alpha_c_x  = 0;
    alpha_c_y  = 0;
    alpha_c_z  = 0;
    alpha_c_dp = 0;

    gamma_c_x  = 0;
    gamma_c_y  = 0;
    gamma_c_z  = 0;
    gamma_c_dp = 0;

    couuang = 0;

    t61 = 0;
    t62 = 0;
    t63 = 0;
    t64 = 0;

    t11 = 0;
    t12 = 0;
    t13 = 0;
    t14 = 0;
}

void OpticsRootOutput::WriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alxi,double alxii,double alzi, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double phxpii, double phzpi, double phzpii,  double couuang_in, double t61_in, double t62_in, double t63_in, double t64_in, double t11_in, double t12_in, double t13_in, double t14_in)
{
    mu_x = phi1;
    mu_y = phi2;

    beta_c_x = bexi;
    beta_c_y = bexii;
    beta_c_z = bezi;
    beta_c_dp = bezii;

    alpha_c_x = alxi;
    alpha_c_y = alxii;
    alpha_c_z = alzi;
    alpha_c_dp = alzii;

    gamma_c_x = gaxi;
    gamma_c_y = gaxii;
    gamma_c_z = gazi;
    gamma_c_dp = gazii;

    phi_x = phxi;
    phi_y = phxii;
    phi_z = phzi;
    phi_dp = phzii;

    phi_xp = phxpi;
    phi_yp = phxpii;
    phi_zp = phzpi;
    phi_dpp = phzpii;

    couuang = couuang_in;

    t11 = t11_in;
    t12 = t12_in;
    t13 = t13_in;
    t14 = t14_in;

    t61 = t61_in;
    t62 = t62_in;
    t63 = t63_in;
    t64 = t64_in;
}

void OpticsRootOutput::WriteLin(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in)
{
    n = n_in;
    std::string tmpname(name_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(name,tmpname.substr(0,c_len).c_str(),49);

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
}

void OpticsRootOutput::Fill()
{
    //Do the write
    OpticsTree->Fill();
}

