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

/*
	if(Track->GetParticleDefinition()->GetAtomicMass() > 1)
	{
		std::cout << "GetCharge() - pre: " << Track->GetDynamicParticle()->GetCharge() << std::endl;
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
//		std::cout << "PARTICLE INTERACTED" << std::endl;
		EventAction->OutputParticle->interacted = 1;
	}
*/

	G4StepStatus Tstatus = Track->GetStep()->GetPostStepPoint()->GetStepStatus();
//    if (Tstatus == fWorldBoundary && Track->GetParticleDefinition() == G4Proton::ProtonDefinition())

    if (Tstatus == fWorldBoundary && Track->GetParticleDefinition()->GetPDGCharge() != 0)
	{
//		std::cout << "AT EXIT PLANE" << std::endl;
		if(Track->GetKineticEnergy() > ReferenceEnergy*RelativeEnergyCut && Track->GetKineticEnergy() > AbsoluteEnergyCut)
		{
			G4Stuff exit_particle;

			exit_particle.x = Track->GetPosition().x()  / CLHEP::m;
			exit_particle.y = Track->GetPosition().y()  / CLHEP::m;
			exit_particle.px = Track->GetMomentum().x() / CLHEP::GeV;
			exit_particle.py = Track->GetMomentum().y() / CLHEP::GeV;
			exit_particle.p  = sqrt(pow(Track->GetMomentum().z()/CLHEP::GeV,2) + pow(Track->GetParticleDefinition()->GetPDGMass()/CLHEP::GeV,2));
			exit_particle.pdgid = Track->GetParticleDefinition()->GetPDGEncoding();
			exit_particle.z = Track->GetParticleDefinition()->GetAtomicNumber();
			exit_particle.a = Track->GetParticleDefinition()->GetAtomicMass();
			exit_particle.m = Track->GetParticleDefinition()->GetPDGMass()/CLHEP::GeV;

			exit_particle.q = Track->GetDynamicParticle()->GetCharge();
			EventAction->AddOutputParticle(exit_particle);
//			std::cout << "TrackEnd: " << Track->GetParticleDefinition()->GetParticleName() << "\t" << Track->GetKineticEnergy() /CLHEP::GeV << "\t" << Track->GetMomentum().z()/CLHEP::GeV << std::endl;
/*
			if(exit_particle.z > 1 && exit_particle.a > 1)
			{
				std::cout << "GetCharge() - post: " << Track->GetDynamicParticle()->GetCharge() << std::endl;
			}
*/
			//double p = Track->GetKineticEnergy();
//			std::cout << "KEEPING: ";// << std::endl;
//			std::cout << x / CLHEP::m << "\t" << px << "\t" << y << "\t" << py << "\t" << p << std::endl;
//			EventAction->AddOutputParticle(x,px,y,py,p);
//      G4double  energy = aStep->GetTrack()->GetKineticEnergy();
//			EventAction->IncrementProtonCount();
		
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

