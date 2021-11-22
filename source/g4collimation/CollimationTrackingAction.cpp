#include "CollimationTrackingAction.h"

#include <iostream>
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

#include "G4Proton.hh"

CollimationTrackingAction::CollimationTrackingAction() : do_debug(false),KeepNeutrals(false),KeepOnlyStable(false),parentID(0),partID_max(0)
{}

void CollimationTrackingAction::PreUserTrackingAction(const G4Track*)
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

//Extraction plane and charge cut
    //if (Tstatus == fWorldBoundary && Track->GetParticleDefinition()->GetPDGCharge() != 0)
    if (Tstatus == fWorldBoundary)
	{
		bool keep_this = false;

//PDG id cuts
		if(keep_ids->empty())
		{
			keep_this = true;
		}
		//Deal with ions
		else if(Track->GetParticleDefinition()->GetPDGEncoding() > 1000000000 && (keep_ids->count(-1) != 0))
		{
			keep_this = true;
		}
		//Non-ions
		else if(keep_ids->count( Track->GetParticleDefinition()->GetPDGEncoding() )  != 0)
		{
			keep_this = true;
		}
		else
		{
			if(do_debug)
			{
				std::cout << "RETURN PID> Not returning particle: " << Track->GetParticleDefinition()->GetParticleName() << "\t" << Track->GetParticleDefinition()->GetPDGEncoding() << std::endl;
			}
			keep_this = false;
		}

//Particle Stability cut
		if(KeepOnlyStable && !Track->GetParticleDefinition()->GetPDGStable())
		{
			keep_this = false;
			if(do_debug)
			{
				std::cout << "RETURN STABLE> Not returning particle: " << Track->GetParticleDefinition()->GetParticleName() << "\t" << Track->GetParticleDefinition()->GetPDGEncoding() << "\t" << !Track->GetParticleDefinition()->GetPDGStable() << std::endl;
			}
		}

//Check if the particle is neutral and if we are keeping them or not
		if(Track->GetParticleDefinition()->GetPDGCharge() != 0 && KeepNeutrals != true)
		{
			keep_this = false;
		}

//Energy cut
		if((Track->GetKineticEnergy() > ReferenceEnergy*RelativeEnergyCut) && (Track->GetKineticEnergy() > AbsoluteEnergyCut) && (keep_this == true))
		{
			CollimationParticle exit_particle;

			exit_particle.x = Track->GetPosition().x()  / CLHEP::m;
			exit_particle.y = Track->GetPosition().y()  / CLHEP::m;
			exit_particle.px = Track->GetMomentum().x() / CLHEP::GeV;
			exit_particle.py = Track->GetMomentum().y() / CLHEP::GeV;
			exit_particle.pz = Track->GetMomentum().z() / CLHEP::GeV;

			double p2 = pow(Track->GetMomentum().x()/CLHEP::GeV,2) + pow(Track->GetMomentum().y()/CLHEP::GeV,2) + pow(Track->GetMomentum().z()/CLHEP::GeV,2); 

			exit_particle.p = sqrt(p2);
			exit_particle.e = sqrt(p2 + pow(Track->GetParticleDefinition()->GetPDGMass()/CLHEP::GeV,2));
			exit_particle.pdgid = Track->GetParticleDefinition()->GetPDGEncoding();
			exit_particle.z = Track->GetParticleDefinition()->GetAtomicNumber();
			exit_particle.a = Track->GetParticleDefinition()->GetAtomicMass();
			exit_particle.m = Track->GetParticleDefinition()->GetPDGMass();

			//Extract the spin vector
			exit_particle.sx = Track->GetPolarization().x();
			exit_particle.sy = Track->GetPolarization().y();
			exit_particle.sz = Track->GetPolarization().z();

			//Extract the time since the initial particle was injected
			exit_particle.t = Track->GetGlobalTime()/CLHEP::second;

			exit_particle.q = Track->GetDynamicParticle()->GetCharge();

			//if the particle is the parent, keep the same parent ID
			//if it is a secondary, use the primary paricle ID as the parent ID
			if(Track->GetParentID() == 0)
			{
				//Do not change anything - keep previous IDs
				exit_particle.parent_id = parentID;
				exit_particle.id        = particleID;
			}
			else
			{
				//Set parent id: use the primary particle ID for this
				exit_particle.parent_id = particleID;

				//Calculate the id of this if it is a new particle.
				partID_max++;
				exit_particle.id = partID_max;
				if(do_debug)
				{
					std::cout << "Making new secondary particle with ID: " << exit_particle.id << " and parent " << exit_particle.parent_id << std::endl;
				}
			}

			//Set the weight of this particle
			exit_particle.weight = Track->GetWeight();

			EventAction->AddOutputParticle(exit_particle);

			if(do_debug)
			{
				std::cout << "RETURN OK> "<< Track->GetParticleDefinition()->GetParticleName() << "\t" << Track->GetParticleDefinition()->GetPDGEncoding() << "\t" << Track->GetKineticEnergy() << std::endl;
			}
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

void CollimationTrackingAction::SetDebug(bool flag)
{
	do_debug = flag;
}

void CollimationTrackingAction::SetParticlesToKeep(std::set<int>* iset)
{
	keep_ids = iset;
}

void CollimationTrackingAction::SetKeepStableParticles(bool flag)
{
	KeepOnlyStable = flag;
}

void CollimationTrackingAction::SetKeepNeutrals(bool flag)
{
	KeepNeutrals = flag;
}

void CollimationTrackingAction::SetParticleID(int32_t in)
{
	particleID = in;
}

int32_t CollimationTrackingAction::GetParticleID()
{
	return particleID;
}

void CollimationTrackingAction::SetParentID(int32_t in)
{
	parentID = in;
}

int32_t CollimationTrackingAction::GetParentID()
{
	return parentID;
}

void CollimationTrackingAction::SetMaximumParticleID(int32_t in)
{
	partID_max = in;
}

int32_t CollimationTrackingAction::GetMaximumParticleID()
{
	return partID_max;
}

