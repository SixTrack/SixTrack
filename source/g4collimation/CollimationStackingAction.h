#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"

#include <fstream>

class CollimationStackingAction : public G4UserStackingAction
{
public:
	CollimationStackingAction() : do_debug(false) {};
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);

	void SetReferenceEnergy(double e0);
	void SetRelativeEnergyCut(double cut);
	void SetAbsoluteEnergyCut(double cut);
	void SetRigidityCut(double cut);
	void SetDebug(bool flag);

private:

	double ReferenceEnergy;
	double RelativeEnergyCut;
	double AbsoluteEnergyCut;
	double RigidityCut;

	bool do_debug;
};

#endif

