#include "CollimationJawSD.h"
#include "CollimationJawHit.h"
#include "G4SDManager.hh"

CollimationJawSD::CollimationJawSD(const G4String& name, const G4String& hitsCollectionName) : G4VSensitiveDetector(name), JawHitsCollection(nullptr), ThisName(name)
{
	collectionName.insert(hitsCollectionName);
}

CollimationJawSD::~CollimationJawSD() {}

G4bool CollimationJawSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) 
{
	double edep = aStep->GetTotalEnergyDeposit();

	if(edep == 0)
	{
		return false;
	}

	CollimationJawHit* ThisHit = new CollimationJawHit();

//	ThisHit->SetTrackID(aStep->GetTrack()->GetTrackID());
//	ThisHit->SetPDGid(aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
	ThisHit->SetEdep(edep);
	ThisHit->SetStartPosition(aStep->GetPreStepPoint()->GetPosition());
	ThisHit->SetEndPosition(aStep->GetPostStepPoint()->GetPosition());
	JawHitsCollection->insert(ThisHit);

	return true;
}

void CollimationJawSD::Initialize(G4HCofThisEvent* hc)
{
	//Set up, etc
	JawHitsCollection = new CollimationJawHitsCollection(SensitiveDetectorName, collectionName[0]);

	int id = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	hc->AddHitsCollection(id, JawHitsCollection);
}

void CollimationJawSD::EndOfEvent(G4HCofThisEvent*)
{
/*
	//Print/historgram
	if(JawHitsCollection->entries() != 0)
	{
		std::cout << "GEANT4> End of event " << ThisName << ": - JawHitsCollection.entries() " <<  JawHitsCollection->entries() << std::endl; 
	}
*/
}

