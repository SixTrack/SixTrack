void canonical2emittance_(double cancord[6], double emittance[3]);
void dist2sixcoord_();
void action2canonical_(double acangl[6], double cancord[6], double acoord[6]);
void action2sixinternal_(double tc[6], double *results, double *normalized);
int checkdist();
void createrandomdist_();
int particle_within_limits_physical(double *physical);
int particle_within_limits_normalized(double *normalized);