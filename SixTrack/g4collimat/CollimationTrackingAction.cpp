#include "CollimationTrackingAction.h"

#include <iostream>
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

#include "G4Proton.hh"

CollimationTrackingAction::CollimationTrackingAction()
{}

void CollimationTrackingAction::PreUserTrackingAction(const G4Track* Track)
{
/*
	const G4ParticleDefinition* particle = Track->GetParticleDefinition();
	if(particle == G4Proton::ProtonDefinition())
	{
		const G4VProcess* pr = Track->GetCreatorProcess();
		if(pr)
		{
			if(pr->GetProcessName() == "protonInelastic")
			{
				double KE = Track->GetKineticEnergy()/CLHEP::GeV;
				if(KE > 1000)
				{
					std::cout << "PROTON INELASTIC: " << KE << std::endl;
					std::cout << Track->GetTrackID() << "\t" << Track->GetParentID() << std::endl;
					//const G4ThreeVector& pos = Track->GetPreStepPoint();
					//std::cout << pos / CLHEP::m << std::endl;
				}
			}
		}
	}
	else
	{
		//Track->SetTrackStatus(fStopAndKill);
	}
*/

//G4StepStatus Tstatus = Track->GetStep()->GetPreStepPoint()->GetStepStatus();
//if(Tstatus == fGeomBoundary)
//Record initial impact point
}

void CollimationTrackingAction::PostUserTrackingAction(const G4Track* Track)
{
	//Check if we have impacted
	G4StepPoint* point_in = Track->GetStep()->GetPreStepPoint();
	G4TouchableHandle touch_in = point_in->GetTouchableHandle();
	G4VPhysicalVolume* volume = touch_in->GetVolume();
	if(volume->GetName().substr(0,3) == "jaw")
	{
//		std::cout << "PARTICLE INTERACTED" << std::endl;
		EventAction->OutputParticle->interacted = 1;
	}

	G4StepStatus Tstatus = Track->GetStep()->GetPostStepPoint()->GetStepStatus();
    if (Tstatus == fWorldBoundary && Track->GetParticleDefinition() == G4Proton::ProtonDefinition())
	{
//		std::cout << "AT EXIT PLANE" << std::endl;
		if(Track->GetKineticEnergy() > ReferenceEnergy*EnergyCut)
		{
			double x = Track->GetPosition().x();
			double y = Track->GetPosition().y();
			double px = Track->GetMomentum().x();
			double py = Track->GetMomentum().y();
			double p = Track->GetMomentum().z();
//			std::cout << "KEEPING: ";// << std::endl;
//			std::cout << x / CLHEP::m << "\t" << px << "\t" << y << "\t" << py << "\t" << p << std::endl;
			EventAction->SetOutputParticle(x,px,y,py,p);
//      G4double  energy = aStep->GetTrack()->GetKineticEnergy();
			EventAction->IncrementProtonCount();
		}
	}
	
}

void CollimationTrackingAction::SetEventAction(CollimationEventAction* ev)
{
	EventAction = ev;
}

void CollimationTrackingAction::SetReferenceEnergy(double e0)
{
	ReferenceEnergy = e0;
}

void CollimationTrackingAction::SetEnergyCut(double cut)
{
	EnergyCut = cut;
}
