#include <cmath>
#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include "CollimationEventAction.h"
#include "CollimationJawHit.h"

#ifdef USE_ROOT_FLAG
#include "RootEnergyDeposition.h"
#endif

CollimationEventAction::CollimationEventAction(bool edep) : G4UserEventAction(), DoEnergyDeposition(edep), ThisCollimatorName(""), ThisCollimatorLength(0.0), ThisCollimatorHalfGap(0.0)
{

}

void CollimationEventAction::BeginOfEventAction(const G4Event* event)
{
//	ThisEvent = event;
}

void CollimationEventAction::EndOfEventAction(const G4Event* ThisEvent)
{
#ifdef USE_ROOT_FLAG
	if(DoEnergyDeposition)
	{
		G4HCofThisEvent* EventHitCollection = ThisEvent->GetHCofThisEvent();
		G4SDManager* sdm = G4SDManager::GetSDMpointer();

		CollimationJawHitsCollection* LeftHits;
		CollimationJawHitsCollection* RightHits;

		int LeftID = sdm->GetCollectionID(ThisCollimatorName + "_LeftJaw");;
		int RightID = sdm->GetCollectionID(ThisCollimatorName + "_RightJaw");

		if(LeftID >= 0)
		{
			LeftHits = (CollimationJawHitsCollection*) EventHitCollection->GetHC(LeftID);
		}
		if(RightID >= 0)
		{
			RightHits = (CollimationJawHitsCollection*) EventHitCollection->GetHC(RightID);
		}

		if(LeftID < 0 || RightID < 0)
		{
			std::cerr << "GEANT4> ERROR: Could not find SD entries for " <<  ThisCollimatorName << std::endl;
			exit(EXIT_FAILURE);
		}

		if(LeftHits && RightHits)
		{
			if(LeftHits->entries() != 0 || RightHits->entries() != 0)
			{
				//Find the histograms for this collimator
				RootOutput->SetCollimator(ThisCollimatorName, ThisCollimatorLength, ThisCollimatorHalfGap);
				//Tell it about the hit collection
				RootOutput->Process(LeftHits);
				RootOutput->Process(RightHits);
			}
		}
	} //End DoEnergyDeposition
#endif
}

void CollimationEventAction::SetOutputVector(std::vector<CollimationParticle>* out)
{
	output_particles = out;
}

void CollimationEventAction::AddOutputParticle(CollimationParticle aParticle)
{
	output_particles->push_back(aParticle);
}

#ifdef USE_ROOT_FLAG
void CollimationEventAction::SetRootOutput(RootEnergyDeposition* root)
{
	RootOutput = root;
}
#endif


void CollimationEventAction::SetCollimator(std::string CollimatorName, double length, double gap)
{
	ThisCollimatorName = CollimatorName;
	ThisCollimatorLength = length;
	ThisCollimatorHalfGap = gap;
}
