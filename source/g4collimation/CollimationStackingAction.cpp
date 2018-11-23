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

	//Always kill things that go below the energy cut
	if(aTrack->GetKineticEnergy() < AbsoluteEnergyCut)
	{
		if(do_debug)
		{
			std::cout << "ABSENERGYCUT> Killing particle with Kinetic energy " << aTrack->GetKineticEnergy() << " < " << AbsoluteEnergyCut << " MeV" << std::endl;
		}
		return fKill;
	}
	else
	{
		return fUrgent;
	}
}

void CollimationStackingAction::SetReferenceEnergy(double e0)
{
	ReferenceEnergy = e0;
}

void CollimationStackingAction::SetAbsoluteEnergyCut(double cut)
{
	AbsoluteEnergyCut = cut;
}

void CollimationStackingAction::SetRigidityCut(double cut)
{
	RigidityCut = cut;
}

void CollimationStackingAction::SetRelativeEnergyCut(double cut)
{
	RelativeEnergyCut = cut;
}

void CollimationStackingAction::SetDebug(bool flag)
{
	do_debug = flag;
}

