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

CollimationEventAction::CollimationEventAction() : InputParticle(new TempParticle), OutputParticle(new TempParticle)
{}

void CollimationEventAction::BeginOfEventAction(const G4Event* event)
{
	ProtonCount = 0;
	ThisEvent = event;
}

void CollimationEventAction::EndOfEventAction(const G4Event*)
{
//	if(ProtonCount != 0)
//	std::cout << "Exiting were " << ProtonCount << " protons" << std::endl;
}


void CollimationEventAction::IncrementProtonCount()
{
	ProtonCount++;
}

void CollimationEventAction::SetInputParticle(double x, double px, double y, double py, double p)
{
	InputParticle->x = x;
	InputParticle->xp = px;
	InputParticle->y = y;
	InputParticle->yp = py;
	InputParticle->p = p;
	size_t oldprec = std::cout.precision(15);
	//std::cout << "IN:  " << x << "\t" << px << "\t" << y << "\t" << py << "\t" << p << std::endl;
	std::cout.precision(oldprec);
}

void CollimationEventAction::SetOutputParticle(double x, double px, double y, double py, double p)
{
	OutputParticle->x = x;
	OutputParticle->xp = px;
	OutputParticle->y = y;
	OutputParticle->yp = py;

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("proton");
	const G4double mp = particle->GetPDGMass();
	double pp= sqrt(p*p + mp*mp);

	OutputParticle->p = pp;
	size_t oldprec = std::cout.precision(15);
	//std::cout << "OUT: " << x << "\t" << px << "\t" << y << "\t" << py << "\t" << pp << std::endl;
	std::cout.precision(oldprec);
}

void CollimationEventAction::PostProcessEvent(double* x, double* y, double* xp, double* yp, double* p, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract)
{
/*
	double x;
	double xp;
	double y;
	double yp;
	double dp;
	double z_inel;
	int interacted;
	int absorbed;
*/
	//Check if we did SD or inelastic + record this track
	//If postprocessing tools need just inelastic, then if we have 0 protons then set fake location for death
	//return true if we want to keep this track
	if(GetProtonCount() == 0)
	{
		//Kill this track, set inelastic flag
		*part_abs = 1; //1 for inelastic, 4 for SD
		*part_hit = 1;
		*part_indiv = 0;
		*part_impact = 100;
		*x= 99.99e-3;
		*y= 99.99e-3;
		*xp= 0.0;
		*yp= 0.0;
		*p= 1;
	}
	else if(GetProtonCount() == 1)
	{
		*part_abs = 0; //1 for inelastic, 4 for SD
		*part_indiv = 0;
		*part_impact = 100;
		*part_linteract = 0.1;
		*x= OutputParticle->x / CLHEP::m;
		*y= OutputParticle->y / CLHEP::m;
		//*xp= tan(OutputParticle->xp / OutputParticle->p);
		//*yp= tan(OutputParticle->yp / OutputParticle->p);
		*xp= (OutputParticle->xp / OutputParticle->p);
		*yp= (OutputParticle->yp / OutputParticle->p);
		*p = OutputParticle->p / CLHEP::GeV;

		//std::cout << "SURVIVED: " << OutputParticle->x << "\t" << OutputParticle->y << std::endl;
		if(OutputParticle->interacted == 1)
		{
			*part_hit = 1;
		}
	}
	else
	{
		std::cerr << "ERROR: More than 1 proton above the energy cut exited the collimator, this means something has gone wrong in geant4" << std::endl;
		exit(EXIT_FAILURE);
	}

	//Finally reset
	if(InputParticle)
	{
		delete InputParticle;
	}
	InputParticle = new TempParticle;
	InputParticle->x = 0;
	InputParticle->xp = 0;
	InputParticle->y = 0;
	InputParticle->yp = 0;
	InputParticle->p = 0;
	InputParticle->interacted = 0;

	if(OutputParticle)
	{
		delete OutputParticle;
	}
	OutputParticle = new TempParticle;
	OutputParticle->x = 0;
	OutputParticle->xp = 0;
	OutputParticle->y = 0;
	OutputParticle->yp = 0;
	OutputParticle->interacted = 0;
}

unsigned int CollimationEventAction::GetProtonCount() const
{
	return ProtonCount;
}
