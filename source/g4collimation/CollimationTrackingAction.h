#ifndef CollimationTrackingAction_h
#define CollimationTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "CollimationEventAction.h"

#include "G4Track.hh"

#include <set>

class CollimationTrackingAction : public G4UserTrackingAction
{
public:

	CollimationTrackingAction();
	void PreUserTrackingAction(const G4Track*);
	void PostUserTrackingAction(const G4Track*);

	void SetEventAction(CollimationEventAction* ev);

	void SetDebug(bool flag);
	void SetParticlesToKeep(std::set<int>*);
	void SetKeepStableParticles(bool);
	void SetKeepNeutrals(bool);

	bool do_debug;
	bool KeepNeutrals;
	bool KeepOnlyStable;

	CollimationEventAction* EventAction;
	void SetReferenceEnergy(double e0);
	void SetRelativeEnergyCut(double cut);
	void SetAbsoluteEnergyCut(double cut);
	void SetRigidityCut(double cut);

	void SetParticleID(int32_t);
	int32_t GetParticleID();

	void SetParentID(int32_t);
	int32_t GetParentID();

	void SetMaximumParticleID(int32_t);
	int32_t GetMaximumParticleID();

	//This particle's ID
	int32_t particleID;

	//This particle's parent particle ID
	int32_t parentID;

	//Maximum partID
	int32_t partID_max;

	double ReferenceEnergy;
	double RelativeEnergyCut;
	double AbsoluteEnergyCut;
	double RigidityCut;

	std::set<int>* keep_ids;

};

#endif

