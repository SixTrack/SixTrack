struct distparam
{
	struct parameters** coord;
    int disttype;
	struct refparam* ref;
	double **tas;
	double **invtas;
	double *closedorbit;
	int coordtype; // This tells which type of coordinates the input is given.  // 1-Normalized
	struct coordinates** incoord;
	struct coordinates** outcoord;
	int totincoord;
	int totoutcoord;
	int isDistrcalculated;
	struct appliedcut* cuts2apply;
};
struct refparam{
	int charge0;
	int z0;
	int a0;
	double pc0;
	double e0;
	double mass0;
	double beta0;
	int en_like;
	int time_like;
	int ang_like;
	struct emittances* emitt;

};
struct coordinates
{
	double *action;
	double *normalized;
	double *physical;
	double *nonstandard;
	double mass;
	int charge;
	int a;
	int z;
	int typeused; // This says which is used or a mixture of it..
	// 0-action, 1-normalized, 2-physical, 3-actionangle + deltap and deltas  
};

struct parameters
{
  double start;
  double stop;   
  int length;
  int type; //This gives the type of distribution, constant, linear, gaussian, 
  double * values;
  int coordtype; 
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

void initializedistribution_(int numberOfDist);
void setdistribution(int ndist);
void setphysicalcut(int variable, double min, double max);
void setnormalizedcut(int variable, double min, double max);



