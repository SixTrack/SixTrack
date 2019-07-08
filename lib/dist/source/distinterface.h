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



