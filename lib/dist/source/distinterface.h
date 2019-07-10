struct cut{
	int isset;
	double min;
	double max;
};

struct distparam* dist;
struct distparam* diststart;
int dim;
int distn;
int writefile_f(const char*  filename_in, int strlen);
int readfile_f(const char*  filename_in, int strlen);
void initializedistribution(int numberOfDist);
void setdistribution(int ndist);
void setphysicalcut(int variable, double min, double max);
void setnormalizedcut(int variable, double min, double max);

void sete0andmass0(double energy0, double mass0);
void setemitt12(double e1, double e2);
void setemitt3(double e3);
void settasmatrix(double * tas);
void setnormalizedcoords(double *xn, double *xnp, double *yn, double *ynp, double *zn, double *znp, int *totparticles);
void get6trackcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles);