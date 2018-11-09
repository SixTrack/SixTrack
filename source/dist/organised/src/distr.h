struct distparam
{
	struct parameters** coord;
	struct emittances* emitt;
	double mass;	
	double momentum;
	double **tas;
	double **invtas;
	int coordtype; // This tells which type of coordinates the input is given.  // 1-Normalized 
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
};


void action2canonical_(double tc[6], double cancord[6]);
int getnumberdist_();
void setemittance12_(double *e1, double *e2);
void initializedistribution_(int *numberOfDist, int *dimension);
void printdistsettings_(int *ndist);
void settasmatrix_(double tas[6][6]);
void dist2sixcoord_(double **results);
void setmassmom_(double *mass, double *momentum);
void setparameter_(int *index,  double *start, double *stop, int *length, int *type);
void deltap(double *dp);
double convertdp2emittance(double dp, double **invtas);
void setemittance3_(double *e3);
void createTasWithNoCoupling(double betax, double alfax, double betay, double alfay, double tas[6][6]);
void action2sixinternal_(double tc[6], double results[6]);
int checkdist();