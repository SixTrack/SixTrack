#ifndef Optics_root_h
#define Optics_root_h

extern "C" void OpticsRootInit();
extern "C" void OpticsRootWriteLin(int n_in, char* name_in, int c_len, double s_in, double x_in, double xp_in, double y_in, double yp_in, double beta_x_in, double beta_y_in, double alpha_x_in, double alpha_y_in, double dispersion_x_in, double dispersion_y_in, double dispersion_xp_in, double dispersion_yp_in);
extern "C" void OpticsRootWriteCpl(double phi1, double phi2, double bexi, double bexii, double bezi, double bezii, double alxi,double alxii,double alzi, double alzii, double gaxi, double gaxii, double gazi, double gazii, double phxi, double phxii, double phzi, double phzii, double phxpi, double phxpii, double phzpi, double phzpii,  double couuang_in, double t61_in, double t62_in, double t63_in, double t64_in, double t11_in, double t12_in, double t13_in, double t14_in);
extern "C" void OpticsRootWrite();

#endif

