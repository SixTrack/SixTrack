#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"

#include <fstream>

class CollimationStackingAction : public G4UserStackingAction
{
public:
	CollimationStackingAction() {};
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
};

#endif

