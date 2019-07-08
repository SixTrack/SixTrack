#ifndef CollimationEventAction_h
#define CollimationEventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"

#include "CollimationStorage.h"

class CollimationEventAction : public G4UserEventAction
{
public:

	CollimationEventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void SetOutputVector(std::vector<CollimationParticle>*);
	void AddOutputParticle(CollimationParticle);

private:
	unsigned int ProtonCount;
	const G4Event* ThisEvent;
	std::vector<CollimationParticle>* output_particles;
};

#endif

