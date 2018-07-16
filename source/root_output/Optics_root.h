#ifndef Optics_root_h
#define Optics_root_h

//Interface calls
extern "C" void OpticsRootInit();
extern "C" void OpticsRootWriteLin(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in);
extern "C" void OpticsRootWriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alxi,double alxii,double alzi, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double phxpii, double phzpi, double phzpii,  double couuang_in, double t61_in, double t62_in, double t63_in, double t64_in, double t11_in, double t12_in, double t13_in, double t14_in);
extern "C" void OpticsRootWrite();

class OpticsRootOutput
{
public:

//Constructor
    OpticsRootOutput();

//Fill optics parameters
    void WriteLin(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in);
    void WriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alxi,double alxii,double alzi, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double phxpii, double phzpi, double phzpii,  double couuang_in, double t61_in, double t62_in, double t63_in, double t64_in, double t11_in, double t12_in, double t13_in, double t14_in);

//Do the actual Fill
    void Fill();

private:

    //Storage tree
    TTree *OpticsTree;


//extern "C" void OpticsRootWriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double couuang, double t61, double t62, double t63, double t64, double t11, double t12, double t13, double t14)
//Must fill the trees from different functions.

/**
subroutine writelin(nr,typ,tl,p1,t,ixwl,isBLOC)

and

subroutine cpltwis(typ,t,etl,phi)
!-----------------------------------------------------------------------
!  CALCULATES COUPLED TWISS PARAMETERS AROUND THE RING AND ALSO THE
!  ANGLE OF THE MAJOR AXIS OF A ELLIPSE IN THE X-Y PROJECTION WITH
!  THE X-AXIS. THE 4-D ELLIPSOID IS GIVEN BY THE BOUNDARY OF A
!  DISTRIBUTION OF PARTICLES WITH MAXIMUM EMITANCE OF MODE I AND II,
!  EUI AND EUII RESPECTIVELY.
!  BINARY PRINT ON FILE 11 OF 22 VALUES :
!  POSITION [M],
!  BET(1-4), ALF(1-4), GAM(1-4), COOR-PHI(1-4), COOR-PRIME-PHI(1-4),
!  COUUANGL
!-----------------------------------------------------------------------
*/

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
Double_t gamma_x;
Double_t gamma_y;


Double_t dispersion_x;
Double_t dispersion_y;

Double_t dispersion_xp;
Double_t dispersion_yp;

Double_t phi_x;
Double_t phi_y;
Double_t phi_z;
Double_t phi_dp;

Double_t phi_xp;
Double_t phi_yp;
Double_t phi_zp;
Double_t phi_dpp;

Double_t mu_x;
Double_t mu_y;

Double_t beta_c_x;
Double_t beta_c_y;
Double_t beta_c_z;
Double_t beta_c_dp;

Double_t alpha_c_x;
Double_t alpha_c_y;
Double_t alpha_c_z;
Double_t alpha_c_dp;

Double_t gamma_c_x;
Double_t gamma_c_y;
Double_t gamma_c_z;
Double_t gamma_c_dp;

Double_t couuang;

Double_t t61;
Double_t t62;
Double_t t63;
Double_t t64;

Double_t t11;
Double_t t12;
Double_t t13;
Double_t t14;

Char_t name[49];
};

#endif

