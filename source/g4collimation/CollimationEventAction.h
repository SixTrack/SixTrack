#ifndef CollimationEventAction_h
#define CollimationEventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"

#ifdef USE_ROOT_FLAG
#include "RootEnergyDeposition.h"
#endif

#include "CollimationStorage.h"

class CollimationEventAction : public G4UserEventAction
{
public:

	CollimationEventAction(bool);
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void SetOutputVector(std::vector<CollimationParticle>*);
	void AddOutputParticle(CollimationParticle);

	void SetCollimator(std::string, double, double);

#ifdef USE_ROOT_FLAG
	void SetRootOutput(RootEnergyDeposition* root);
#endif

private:

	bool DoEnergyDeposition;

	std::string ThisCollimatorName;
	double ThisCollimatorLength;
	double ThisCollimatorHalfGap;

	std::vector<CollimationParticle>* output_particles;

#ifdef USE_ROOT_FLAG
	RootEnergyDeposition* RootOutput;
#endif

};

#endif

