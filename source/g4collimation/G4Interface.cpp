#include <iostream>

#include "G4RunManager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"

#include "CollimationGeometry.h"
#include "CollimationParticleGun.h"
#include "CollimationStackingAction.h"
#include "CollimationTrackingAction.h"

#include <string>
#include <set>

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "CollimationStorage.h"

CollimationParticleGun* part;
G4RunManager* runManager;
//CollimationGeometry* CollimatorJaw;
FTFP_BERT* physlist_FTFP;
QGSP_BERT* physlist_QGSP;
CollimationStackingAction* stack;
CollimationEventAction* event;
CollimationTrackingAction* tracking;
CollimationGeometry* geometry;


std::string CleanFortranString(char* str, size_t count);

std::vector<G4Stuff> input_particles;
std::vector<G4Stuff> output_particles;

std::set<int> keep_ids;

/**
Geant4 needs the user to define several user classes.
These include:
1: The geometry (our collimator jaws). See the local "G4Collimator" class
2: The Physics to use.
3: A particle source.
*/
extern "C" void g4_collimation_init_(double* ReferenceE, int* seed, double* recut, double* aecut, double* rcut, int* PhysicsSelect, bool* g4_debug, bool* g4_keep_stable)
{
	std::cout << "Using seed " << *seed << " in geant4 C++" << std::endl;
	std::cout << "The reference energy is " << *ReferenceE / CLHEP::GeV<< " and the relative energy cut will be at "<< (*ReferenceE * *recut ) / CLHEP::GeV << " GeV!" << std::endl;
	std::cout << "The reference energy is " << *ReferenceE / CLHEP::GeV<< " and the absolute energy cut will be at "<< *aecut << " GeV!" << std::endl;
	CLHEP::HepRandom::setTheSeed(*seed);

	//Construct the run manager
	runManager = new G4RunManager();

	//Physics list
	G4int verbose = 0;
	if(*PhysicsSelect == 0)
	{
		physlist_FTFP = new FTFP_BERT(verbose);
		//Add the physics
		runManager->SetUserInitialization(physlist_FTFP);
	}
	else if(*PhysicsSelect == 1)
	{
		physlist_QGSP = new QGSP_BERT(verbose);
		//Add the physics
		runManager->SetUserInitialization(physlist_QGSP);
	}
	else
	{
		std::cerr << "ERROR: Bad value for Geant4 physics list selection: Exiting" << std::endl;
		exit(EXIT_FAILURE);
	}

	//Construct our collimator jaw geometry
	geometry = new CollimationGeometry();
	geometry->SetDebug(*g4_debug);

	//Make the particle gun
	part = new CollimationParticleGun();

	//This is in MeV in both sixtrack (e0) and in geant4.
	part->SetReferenceEnergy(*ReferenceE);
	part->SetDebug(*g4_debug);

	event = new CollimationEventAction();
	event->SetOutputVector(&output_particles);

	tracking = new CollimationTrackingAction();
	tracking->SetEventAction(event);
	tracking->SetReferenceEnergy(*ReferenceE);
	tracking->SetRelativeEnergyCut(*recut);
	tracking->SetAbsoluteEnergyCut(*aecut * CLHEP::GeV);
	tracking->SetRigidityCut(*rcut);

	if(*g4_debug)
	{
		std::cout << "GEANT4> Keeping: " << keep_ids.size() << " types of particles" << std::endl;
		std::set<int>::const_iterator it = keep_ids.begin();
		while(it != keep_ids.end())
		{
			std::cout << "GEANT4 RETURNING> " << *it << std::endl;
			it++;
		}
	}
	tracking->SetParticlesToKeep(&keep_ids);

	tracking->SetKeepStableParticles(g4_keep_stable);
	tracking->SetDebug(*g4_debug);


	stack = new CollimationStackingAction();
	stack->SetDebug(*g4_debug);
	stack->SetReferenceEnergy(*ReferenceE);
	stack->SetRelativeEnergyCut(*recut);
	stack->SetAbsoluteEnergyCut(*aecut * CLHEP::GeV);
	stack->SetRigidityCut(*rcut);

	//Add our geometry class
	runManager->SetUserInitialization(geometry);

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

	std::string CollimatorName = CleanFortranString(name, 48);
	std::string CollimatorMaterialName = CleanFortranString(material, 4);
	std::cout << "Adding \"" << CollimatorName << "\" with material \"" << CollimatorMaterialName << "\" and rotation \"" << *rotation << "\" and offset \"" << *offset << "\" and length \"";
	std::cout << *length << "\"" << std::endl;

	G4double length_in = *length *CLHEP::m;
	G4double aperture_in = *aperture *CLHEP::m;
	G4double rotation_in = *rotation *CLHEP::rad;
	G4double offset_in = *offset *CLHEP::m;

	geometry->AddCollimator(CollimatorName, length_in, aperture_in, rotation_in, offset_in, CollimatorMaterialName);

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
	std::string CollimatorName = CleanFortranString(name, 48);
	geometry->SetCollimator(CollimatorName);

	G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

	runManager->ReinitializeGeometry();
}

extern "C" void g4_add_particle_(double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, int16_t* nzz, int16_t* naa, int16_t* nqq, double* mass)
{
//WARNING: at this stage in SixTrack the units have been converted to GeV, m, and rad!
//The particle energy input is the TOTAL energy
//mass (i.e. nucm) is already in MeV!

	double x_in = (*x) * CLHEP::m;
	double y_in = (*y) * CLHEP::m;

//We want px and py, not the angle!
	double e_in = (*e) * CLHEP::GeV;
	double p_in = sqrt((e_in*e_in) - (*mass * *mass));

	double px_in = p_in * (*xp);
	double py_in = p_in * (*yp);

//	double px_in = p_in * (*xp) * CLHEP::GeV;
//	double py_in = p_in * (*yp) * CLHEP::GeV;

	G4Stuff in_particle;
	in_particle.x = x_in;
	in_particle.y = y_in;

	in_particle.px = px_in;
	in_particle.py = py_in;
	in_particle.p = p_in;

	in_particle.e = e_in;
	in_particle.pdgid = *pdgid;
	in_particle.z = *nzz;
	in_particle.a = *naa;
	in_particle.q = *nqq;
	in_particle.id = input_particles.size();

	input_particles.push_back(in_particle);
}

extern "C" void g4_collimate_()
{
	output_particles.clear();
	//Update the gun with this particle's details
	for(size_t n=0; n < input_particles.size(); n++)
	{
		part->SetParticleDetails(input_particles.at(n).x, input_particles.at(n).y, input_particles.at(n).px, input_particles.at(n).py, input_particles.at(n).e, input_particles.at(n).p, input_particles.at(n).pdgid, input_particles.at(n).q);

		//Run!
		runManager->BeamOn(1);
	}
	input_particles.clear();
}

/**
* Here we put the particles back into sixtrack and set any flags if needed
*/
extern "C" void g4_collimate_return_(int* j, double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, double* m, int16_t* z, int16_t* a, int16_t* q, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract)
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
//event->PostProcessEvent(x, y, xp, yp, p, part_hit, part_abs, part_impact, part_indiv, part_linteract);
//CollimatorKeyMap.clear();
/*
*x= OutputParticle->x / CLHEP::m;
*y= OutputParticle->y / CLHEP::m;
*xp= (OutputParticle->xp / OutputParticle->p);
*yp= (OutputParticle->yp / OutputParticle->p);
*p = OutputParticle->p / CLHEP::GeV;
*/
*x  = output_particles.at(*j).x;
*y  = output_particles.at(*j).y;
*xp = output_particles.at(*j).px / output_particles.at(*j).p;
*yp = output_particles.at(*j).py / output_particles.at(*j).p;
*e  = output_particles.at(*j).e;
*pdgid  = output_particles.at(*j).pdgid;
*z = output_particles.at(*j).z;
*a  = output_particles.at(*j).a;
*q  = output_particles.at(*j).q;
*m  = output_particles.at(*j).m;
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

extern "C" void g4_get_particle_count_(int* g4_npart)
{
	*g4_npart = output_particles.size(); 
}

extern "C" void g4_collimation_clear_()
{
	input_particles.clear();
	output_particles.clear();
}

extern "C" void g4_keep_id_(int* id)
{
	keep_ids.insert(*id);
}
