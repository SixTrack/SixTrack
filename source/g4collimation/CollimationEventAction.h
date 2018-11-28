#ifndef CollimationEventAction_h
#define CollimationEventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"

#include "Storage.h"

class CollimationEventAction : public G4UserEventAction
{
public:

	CollimationEventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void SetOutputVector(std::vector<G4Stuff>*);
	void AddOutputParticle(G4Stuff);

private:
	unsigned int ProtonCount;
	const G4Event* ThisEvent;
	std::vector<G4Stuff>* output_particles;
};

#endif

