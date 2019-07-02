#define MAX_COLUMNS 100
#define MAX_LENGTH 100
#define MAX_ROWS 100000
void canonical2six(double *canonical, double beta0, double pc0, double mass0, double mass, double *coord);
double momentum2energy(double momentum, double mass);
int readfile();
int splitline(char* line, char columns[MAX_COLUMNS][MAX_LENGTH],char units[MAX_COLUMNS][MAX_LENGTH] );
void add2table(double ** table, char* line, int linenum );
void allocateincoord(int linecount);
void setphysical(int coordorder, int column, double ** table, double multifactor);
void setaction(int coordorder, int column, double ** table);
void setnormalized(int coordorder, int column, double ** table);
void setnonstandard(int coordorder, int column, double ** table);
void checkifenergyset(int entype);
double energy2momentum(double energy, double mass);
void calculaterefparam();
void convert2standard();
double pt2deltap(double pt, double beta0 );
double psigma2deltap(double psigma, double beta0 );
void canonical2six(double *canonical, double beta0, double pc0, double mass0, double mass, double *coord);
double getEnergyUnit(char *unit);
double getMetricUnit(char *unit);
void issue_error(const char* t1);
void issue_warning(const char* t1);
void issue_info(const char* t1);
void add2internaltab(char* word, char columns[MAX_COLUMNS][MAX_LENGTH], char units[MAX_COLUMNS][MAX_LENGTH], int count) ;
int strcmpnl (const char *s1, const char *s2);
void issue_error2(const char* t1,const char* t2);