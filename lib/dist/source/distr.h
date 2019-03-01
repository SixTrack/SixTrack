struct distparam
{
	struct parameters** coord;
	struct emittances* emitt;
	double mass;	
	double momentum;
	double **tas;
	double **invtas;
	double *closedorbit;
	int coordtype; // This tells which type of coordinates the input is given.  // 1-Normalized
	double **distout;
	int isDistrcalculated;
	int longitunalemittance; // 0 - no longitudnal, 1 - e3, 2 - dp, 3 - deltas
};

struct parameters
{
  double start;
  double stop;   
  int length;
  int type; //This gives the type of distribution, constant, linear, gaussian, 
  double * values; 
};

struct emittances{
	double e1, e2, e3;
	double dp, deltas;
};


void action2canonical_(double tc[6], double cancord[6]);
void setdistribution_(int *ndist);
int getnumberdist_();
void setemittance12_(double *e1, double *e2);
void setemittance3_(double *e3);
void initializedistribution_(int *numberOfDist, int *dimension);
void printdistsettings_(int *ndist);
void addclosedorbit_(double *clo);
void settasmatrix_(double tas[6][6]);
void dist2sixcoord_();
void setmassmom_(double *mass, double *momentum);
void setparameter_(int *index,  double *start, double *stop, int *length, int *type);
void setdeltap_(double *dp);
void convertdp2emittance(double dp);
void setemittance3_(double *e3);
void createTasWithNoCoupling(double betax, double alfax, double betay, double alfay, double tas[6][6]);
void action2sixinternal_(double tc[6], double results[6]);
void calcualteinverse();
int checkdist();
void change_e3_to_dp(double cancord[6],double acoord[6], double acangl[6]);
void change_e3_to_dp_easy(double cancord[6], double acoord[6], double acangl[6]);
double fun (double cancord[6], double acoord[6], double acangl[6], double x);