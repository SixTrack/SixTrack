#ifndef CollimationEventAction_h
#define CollimationEventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"

#include "Storage.h"

class CollimationEventAction : public G4UserEventAction
{
public:

	CollimationEventAction();
	//~CollimationEventAction();
	//virtual ~G4UserEventAction() {;}
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void SetOutputVector(std::vector<G4Stuff>*);
	void AddOutputParticle(G4Stuff);

/*
	void IncrementProtonCount();
	unsigned int GetProtonCount() const;
	void SetInputParticle(double, double, double, double, double);
	void SetOutputParticle(double, double, double, double, double);
	void PostProcessEvent(double* x, double* y, double* xp, double* yp, double* dp, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double* part_linteract);

struct TempParticle
{
	double x;
	double xp;
	double y;
	double yp;
	double p;
	double z_inel;
	int interacted;
	int absorbed;
};

TempParticle* InputParticle;
TempParticle* OutputParticle;
*/
private:
	unsigned int ProtonCount;
	const G4Event* ThisEvent;
	std::vector<G4Stuff>* output_particles;
};

#endif

