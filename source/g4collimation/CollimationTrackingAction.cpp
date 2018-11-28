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
}

void CollimationTrackingAction::PostUserTrackingAction(const G4Track* Track)
{
/*
	//Check if we have impacted
	G4StepPoint* point_in = Track->GetStep()->GetPreStepPoint();
	G4TouchableHandle touch_in = point_in->GetTouchableHandle();
	G4VPhysicalVolume* volume = touch_in->GetVolume();
	if(volume->GetName().substr(0,3) == "jaw")
	{
		EventAction->OutputParticle->interacted = 1;
	}
*/

	G4StepStatus Tstatus = Track->GetStep()->GetPostStepPoint()->GetStepStatus();
    if (Tstatus == fWorldBoundary && Track->GetParticleDefinition()->GetPDGCharge() != 0)
	{
		if(Track->GetKineticEnergy() > ReferenceEnergy*RelativeEnergyCut && Track->GetKineticEnergy() > AbsoluteEnergyCut)
		{
			G4Stuff exit_particle;

			exit_particle.x = Track->GetPosition().x()  / CLHEP::m;
			exit_particle.y = Track->GetPosition().y()  / CLHEP::m;
			exit_particle.px = Track->GetMomentum().x() / CLHEP::GeV;
			exit_particle.py = Track->GetMomentum().y() / CLHEP::GeV;

			double p2 = pow(Track->GetMomentum().x()/CLHEP::GeV,2) + pow(Track->GetMomentum().y()/CLHEP::GeV,2) + pow(Track->GetMomentum().z()/CLHEP::GeV,2); 

			exit_particle.p  = sqrt(p2 + pow(Track->GetParticleDefinition()->GetPDGMass()/CLHEP::GeV,2));
			exit_particle.pdgid = Track->GetParticleDefinition()->GetPDGEncoding();
			exit_particle.z = Track->GetParticleDefinition()->GetAtomicNumber();
			exit_particle.a = Track->GetParticleDefinition()->GetAtomicMass();
			exit_particle.m = Track->GetParticleDefinition()->GetPDGMass()/CLHEP::GeV;

			exit_particle.q = Track->GetDynamicParticle()->GetCharge();
			EventAction->AddOutputParticle(exit_particle);
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

void CollimationTrackingAction::SetAbsoluteEnergyCut(double cut)
{
	AbsoluteEnergyCut = cut;
}


void CollimationTrackingAction::SetRigidityCut(double cut)
{
	RigidityCut = cut;
}


void CollimationTrackingAction::SetRelativeEnergyCut(double cut)
{
	RelativeEnergyCut = cut;
}

