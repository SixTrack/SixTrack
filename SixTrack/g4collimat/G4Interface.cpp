#include <iostream>

#include "G4RunManager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"

#include "CollimationGeometry.h"
#include "CollimationParticleGun.h"
//#include "CollimationSteppingAction.h"
#include "CollimationStackingAction.h"
#include "CollimationTrackingAction.h"
#include "CollimationMaterials.h"

#include <map>
#include <string>


CollimationParticleGun* part;
G4RunManager* runManager;
CollimationGeometry* CollimatorJaw;
FTFP_BERT* physlist;
//QGSP_BERT* physlist;
//CollimationSteppingAction* step;
CollimationStackingAction* stack;
CollimationEventAction* event;
CollimationTrackingAction* tracking;
CollimationMaterials* Mmap;

std::map<std::string,CollimationGeometry*> CollimatorKeyMap;

std::string CleanFortranString(char* str, size_t count);


/**
Geant4 needs the user to define several user classes.
These include:
1: The geometry (our collimator jaws). See the local "G4Collimator" class
2: The Physics to use.
3: A particle source.
*/
extern "C" void g4_collimation_init_(double* ReferenceE, int* seed, double* ecut)
{
	std::cout << "Using seed " << *seed << " in geant4 C++" << std::endl;
	std::cout << "The reference energy is " << *ReferenceE / CLHEP::GeV<< " and the cut will be at "<< (*ReferenceE * *ecut ) / CLHEP::GeV << " GeV!" << std::endl;
	CLHEP::HepRandom::setTheSeed(*seed);

	//Physics list
	G4int verbose = 0;
	physlist = new FTFP_BERT(verbose);
//	physlist = new QGSP_BERT(verbose);

	//Construct the run manager
	runManager = new G4RunManager();

	//Add the physics
	runManager->SetUserInitialization(physlist);

	Mmap = new CollimationMaterials();
	G4Material* JawMaterial = Mmap->GetMaterial("Iner");

	//Construct our collimator jaw geometry
	CollimatorJaw = new CollimationGeometry("TempCollimator",10, 0.0,0,0, JawMaterial);

	//Make the particle gun
	part = new CollimationParticleGun();

	//This is in MeV in both sixtrack (e0) and in geant4.
	part->SetReferenceEnergy(*ReferenceE);

	event = new CollimationEventAction();
	tracking = new CollimationTrackingAction();
	tracking->SetEventAction(event);
	tracking->SetEnergyCut(*ecut);
	tracking->SetReferenceEnergy(*ReferenceE);

	stack = new CollimationStackingAction();

	//Add our geometry class
	runManager->SetUserInitialization(CollimatorJaw);

	//Add our particle source class
	runManager->SetUserAction(part);

	//Set up our event action
	runManager->SetUserAction(event);

	//Set up our tracking action
	runManager->SetUserAction(tracking);

	//Set up the stacking action to kill non-protons
	runManager->SetUserAction(stack);

	//Added everything now set up the run manager
	runManager->Initialize();
}

extern "C" void g4_add_collimator_(char* name, char* material, double* length, double* aperture, double* rotation, double* offset)
{
	//NOTE: The closed orbit offset (e.g. at TCTs due to the crossing angle) should be subtracted in the sixtrack part.
	//rcx(j)  = (xv(1,j)-torbx(ie))/1d3
	//rcy(j)  = (xv(2,j)-torby(ie))/1d3
	//Therefore we do not need to take it into account here...

	std::string CollimatorName = CleanFortranString(name, 16);
	std::string CollimatorMaterialName = CleanFortranString(material, 4);
	std::cout << "Adding \"" << CollimatorName << "\" with material \"" << CollimatorMaterialName << "\" and rotation \"" << *rotation << "\"" << std::endl;

	G4Material* JawMaterial = Mmap->GetMaterial(CollimatorMaterialName);

	G4double length_in = *length *CLHEP::m;
	G4double aperture_in = *aperture *CLHEP::m;
	G4double rotation_in = *rotation *CLHEP::rad;
	G4double offset_in = *offset *CLHEP::m;
	CollimationGeometry* NewCollimator = new CollimationGeometry(CollimatorName, length_in, aperture_in, rotation_in, offset_in, JawMaterial);

	std::pair<std::string,CollimationGeometry*> CollimatorKey;
	CollimatorKey = std::make_pair(CollimatorName,NewCollimator);

	std::pair<std::map< std::string,CollimationGeometry* >::const_iterator,bool > check;
	check = CollimatorKeyMap.insert(CollimatorKey);

	if(check.second == false)
	{
		std::cerr << "ERROR: Multiple definitions of collimator: \"" << CollimatorName << "\"" << std::endl;
		exit(EXIT_FAILURE);
	}
}

extern "C" void g4_terminate_()
{
	std::cout << "EXIT GEANT4" << std::endl;

	if(runManager)
	{
		std::cout << "Deleting geant4 run manager" << std::endl;
		delete runManager;
	}
}

//Set up new collimator
//runManager->ReinitializeGeometry();
extern "C" void g4_set_collimator_(char* name)
{
	std::string CollimatorName = CleanFortranString(name, 16);
//	std::cout << "Looking for collimator \"" << CollimatorName << "\"" << std::endl;
	std::map<std::string,CollimationGeometry*>::iterator CollimatorKey_itr = CollimatorKeyMap.find(CollimatorName);
	if(CollimatorKey_itr == CollimatorKeyMap.end())
	{
		std::cerr << "ERROR: Collimator \"" <<  CollimatorName << "\" not found!" << std::endl;
		exit(EXIT_FAILURE);
	}

	runManager->SetUserInitialization(CollimatorKey_itr->second);
	runManager->ReinitializeGeometry();
}

extern "C" void g4_collimate_(double* x, double* y, double* xp, double* yp, double* p)
{
//WARNING: at this stage in SixTrack the units have been converted to GeV, m, and rad!
//The particle energy input is the TOTAL energy
	double x_in = *x * CLHEP::m;
	double y_in = *y * CLHEP::m;

//We want px and py, not the angle!
	double px_in = *p * tan(*xp) * CLHEP::GeV;
	double py_in = *p * tan(*yp) * CLHEP::GeV;
	double p_in = *p * CLHEP::GeV;

//Store the input particle for a check vs the output particle!
	event->SetInputParticle(x_in, y_in, px_in, py_in, p_in);

	//Update the gun with this particle's details
	part->SetParticleDetails(x_in, y_in, px_in, py_in, p_in);

	//Run!
	runManager->BeamOn(1);
}

/**
* Here we put the particles back into sixtrack and set any flags if needed
*/
extern "C" void g4_collimate_return_(double* x, double* y, double* xp, double* yp, double* p, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract)
{
/*
part_hit(j), part_abs(j), part_impact(j), part_indiv(j),
     & part_linteract(j))
!++  Output information:
!++
!++  PART_HIT(MAX_NPART)     Hit flag for last hit (10000*element# + turn#)
!++  PART_ABS(MAX_NPART)     Abs flag (10000*element# + turn#)
!++  PART_IMPACT(MAX_NPART)  Impact parameter (0 for inner face)
!++  PART_INDIV(MAX_NPART)   Divergence of impacting particles
*/
event->PostProcessEvent(x, y, xp, yp, p, part_hit, part_abs, part_impact, part_indiv, part_linteract);
}

std::string CleanFortranString(char* str, size_t count)
{
	//Fortran string nightmares
	std::string Whitespace(" \t\n\r");
	std::string sstring(str,count);

	size_t lpos = sstring.find_last_not_of(Whitespace);
	if(lpos!=std::string::npos)
	{
	    sstring.erase(lpos+1);
	}

	//Fortran string happy place
	return sstring;
}

