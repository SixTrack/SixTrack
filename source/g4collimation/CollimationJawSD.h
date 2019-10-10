#ifndef _CollimationJawSD_h
#define _CollimationJawSD_h

#include <string>

#include "G4VSensitiveDetector.hh"
#include "CollimationJawHit.h"

class CollimationJawSD : public G4VSensitiveDetector
{
public:
//  In the derived class constructor, name(s) of hits collection(s) which
//  // are made by the sensitive detector must be set to "collectionName" string
//  // vector.
	CollimationJawSD(const G4String& name, const G4String& hitsCollectionName);
	~CollimationJawSD();

	G4bool ProcessHits(G4Step*aStep, G4TouchableHistory*ROhist);
	void Initialize(G4HCofThisEvent* hc);
	void EndOfEvent(G4HCofThisEvent* hc);

private:
	CollimationJawHitsCollection* JawHitsCollection;

	std::string ThisName;
};

#endif

