#define MAX_COLUMNS 100
#define MAX_LENGTH 100
#define MAX_ROWS 100000
void canonical2six(double *canonical, double beta0, double pc0, double mass0, double mass, double *coord);
double momentum2energy(double momentum, double mass);
double energy2momentum(double energy, double mass);
double pt2deltap(double pt, double beta0 );
double psigma2deltap(double psigma, double beta0 );
void issue_error(const char* t1);
void issue_error2(const char* t1,const char* t2);
void issue_warning(const char* t1);
void issue_info(const char* t1);
int strcmpnl (const char *s1, const char *s2);
double sigma2zeta(double sigma, double beta0, double beta);
double tau2zeta(double tau, double beta);
