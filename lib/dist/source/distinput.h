struct distparam
{
	struct parameters** coord;
	struct emittances* emitt;
	double mass;
	int charge;
    int massnum;
    int atomnum;
    int totallength;
    int disttype;	
	double momentum;
	double **tas;
	double **invtas;
	double *closedorbit;
	int coordtype; // This tells which type of coordinates the input is given.  // 1-Normalized
	double **distout;
	double **distout_normalized;
	int isDistrcalculated;
	int longitunalemittance; // 0 - no longitudnal, 1 - e3, 2 - dp, 3 - deltas
	struct appliedcut* cuts2apply;
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

struct appliedcut{
	int isset_p;
	int isset_n;
	struct cut** physical;
	struct cut** normalized;
};

struct cut{
	int isset;
	double min;
	double max;
};

struct distparam* dist;
struct distparam* diststart;
int dim;
int distn;



void initializedistribution_(int *numberOfDist);
void setdistribution_(int *ndist);
void setmassmom_(double *mass, double *momentum);
void setemittance12_(double *e1, double *e2);
void setemittance3_(double *e3);
void usedeltap_();


void addclosedorbit_(double *clo);
void settasmatrix_(double tas[6][6]);
void settasmatrixpython(double **tas);
void settasmatrix_element(double element, int i, int j);

void setparameter_(int *index,  double *start, double *stop, int *length, int *type);

//void convertdp2emittance(double dp);

void createtas0coupling_(double betax, double alfax, double betay, double alfay, double dx, double dpx, double dy, double dpy);

//void change_e3_to_dp(double cancord[6],double acoord[6], double acangl[6]);

double optideltas(double cancord[6], double acoord[6], double acangl[6], double x); //check these...
double toactioncord_(double cancord[6], double acoord[6], double acangl[6]);
void setphysicalcut(int variable, double min, double max);
void setnormalizedcut(int variable, double min, double max);
