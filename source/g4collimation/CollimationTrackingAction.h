#ifndef CollimationTrackingAction_h
#define CollimationTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "CollimationEventAction.h"

#include "G4Track.hh"

class CollimationTrackingAction : public G4UserTrackingAction
{
public:

	CollimationTrackingAction();
	void PreUserTrackingAction(const G4Track*);
	void PostUserTrackingAction(const G4Track*);

	void SetEventAction(CollimationEventAction* ev);

	CollimationEventAction* EventAction;
	void SetReferenceEnergy(double e0);
	void SetRelativeEnergyCut(double cut);
	void SetAbsoluteEnergyCut(double cut);
	void SetRigidityCut(double cut);

	double ReferenceEnergy;
	double RelativeEnergyCut;
	double AbsoluteEnergyCut;
	double RigidityCut;
};

#endif

