#include "CollimationStackingAction.h"

#include "G4Proton.hh"
#include "G4Track.hh"

//CollimationStackingAction::CollimationStackingAction() {}

G4ClassificationOfNewTrack CollimationStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
/*
	//See G4 example B3
	if(aTrack->GetParentID() == 0) return fUrgent;

	//Kill all non-protons
	if(aTrack->GetDefinition() != G4Proton::ProtonDefinition())
	{
		//std::cout << "Killing " << aTrack->GetDefinition()->GetParticleName() << std::endl;
		return fKill;
	}
	else
	{
		return fUrgent;
	}
*/
		return fUrgent;
}
