#include "CollimationEventAction.h"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"

#include <cmath>
#include <iostream>

/**
* Adding the collimator geometry lets us know where to extract the beam
* Adding the current event lets us know what the input particle type was, and therefore what type to keep
* This could be protons, Pb ions, etc...
*/

/*
CollimationSteppingAction::CollimationSteppingAction(CollimatorGeometry* col, CollimatorEventAction* evnt)
:collimator(col), eventaction(evt)					 
{}
*/

CollimationEventAction::CollimationEventAction() 
{}

void CollimationEventAction::BeginOfEventAction(const G4Event* event)
{
	ThisEvent = event;
}

void CollimationEventAction::EndOfEventAction(const G4Event*)
{
//	if(ProtonCount != 0)
//	std::cout << "Exiting were " << ProtonCount << " protons" << std::endl;
}

void CollimationEventAction::SetOutputVector(std::vector<CollimationParticle>* out)
{
	output_particles = out;
}

void CollimationEventAction::AddOutputParticle(CollimationParticle aParticle)
{
	output_particles->push_back(aParticle);
}

