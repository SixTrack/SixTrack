
#include <algorithm>
#include <iostream>
#include <set>
#include <string>

#include "CollimationGeometry.h"
#include "CollimationParticleGun.h"
#include "CollimationStackingAction.h"
#include "CollimationStorage.h"
#include "CollimationTrackingAction.h"

#include "CollimationEMD.h"

#ifdef USE_ROOT_FLAG
#include "RootEnergyDeposition.h"
#endif

#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysListFactory.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"

//G4UserRunAction

CollimationParticleGun* part;
G4RunManager* runManager;
CollimationStackingAction* stack;
CollimationEventAction* event;
CollimationTrackingAction* tracking;
CollimationGeometry* geometry;

#ifdef USE_ROOT_FLAG
RootEnergyDeposition* ro;
#endif

std::string CleanFortranString(char* str, size_t count);

std::vector<CollimationParticle> input_particles;
std::vector<CollimationParticle> output_particles;

std::vector<CollimationEnergyDeposition> EnergyDepositionConfiguration;

std::set<int> keep_ids;

int32_t MaximumParticleID;

/**
Geant4 needs the user to define several user classes.
These include:
1: The geometry (our collimator jaws). See the local "G4Collimator" class
2: The Physics to use.
3: A particle source.
*/
extern "C" void g4_collimation_init(double* ReferenceE, int* seed, double* recut, double* aecut, double* rcut, double* rangecut_mm, double* v0, char* PhysicsSelect, bool* g4_debug, bool* g4_keep_stable, bool *DoEnergyDeposition, bool *g4_neutral)
{

	std::cout << "GEANT4> Using seed " << *seed << " in geant4 C++" << std::endl;
	std::cout << "GEANT4> The reference energy is " << *ReferenceE / CLHEP::GeV<< " and the relative energy cut will be at " << (*ReferenceE * *recut ) / CLHEP::GeV << " GeV!" << std::endl;
	std::cout << "GEANT4> The reference energy is " << *ReferenceE / CLHEP::GeV<< " and the absolute energy cut will be at " << *aecut << " GeV!" << std::endl;

	if(*rangecut_mm != 0)
	{
		std::cout << "GEANT4> Using a user-set range cut of " << *rangecut_mm << " mm!" << std::endl;
	}
	else
	{
		std::cout << "GEANT4> Using the default range cuts!" << std::endl;
	}

	std::cout << "GEANT4> Perform energy deposition?: " << *DoEnergyDeposition << std::endl;
	std::cout << "GEANT4> Energy deposition configuration entries: " << EnergyDepositionConfiguration.size()<< std::endl;

	CLHEP::HepRandom::setTheSeed(*seed);

	//Construct the run manager
	runManager = new G4RunManager();

	//Physics list
	G4int verbose = 0;
	std::string PhysicsListName(PhysicsSelect);

	std::cout << "GEANT4> Selected physics list: \"" << PhysicsListName << "\"" << std::endl;

	G4PhysListFactory* PlistFactory = new G4PhysListFactory();
	G4VModularPhysicsList* PhysicsList = PlistFactory->GetReferencePhysList(PhysicsListName);

	if(!PhysicsList)
	{
		std::cerr << "GEANT4> ERROR: Failed to build physics list: \"" << PhysicsSelect << "\"" << std::endl;

		exit(EXIT_FAILURE);
	}

	PhysicsList->SetVerboseLevel(verbose);

	if(*rangecut_mm != 0)
	{
		PhysicsList->SetDefaultCutValue(*rangecut_mm);
	}
	PhysicsList->RegisterPhysics(new EMDissociation());

	runManager->SetUserInitialization(PhysicsList);

	delete PlistFactory;

	//Construct our collimator jaw geometry
	geometry = new CollimationGeometry(*DoEnergyDeposition);
	geometry->SetDebug(*g4_debug);

	//Make the particle gun
	part = new CollimationParticleGun();

	//This is in MeV in both sixtrack (e0) and in geant4.
	part->SetReferenceEnergy(*ReferenceE);
	part->SetDebug(*g4_debug);

	event = new CollimationEventAction(*DoEnergyDeposition);
	event->SetOutputVector(&output_particles);

#ifdef USE_ROOT_FLAG
	if(*DoEnergyDeposition)
	{
		ro = new RootEnergyDeposition(EnergyDepositionConfiguration);
		event->SetRootOutput(ro);

		std::vector<CollimationEnergyDeposition>::const_iterator itr = EnergyDepositionConfiguration.begin();
		while(itr!=EnergyDepositionConfiguration.end())
		{
			std::cout << "GEANT4> Energy deposition: " << itr->name << "\t";
			if(itr->type == 2)
			{
				std::cout << "HIST2D\t";
			}
			else if(itr->type == 3)
			{
				std::cout << "HIST3D\t";

			}

			std::cout << itr->xmax << "\t" << itr->ymax << "\t" << itr->zmax << "\t";
			std::cout << itr->xstep << "\t" << itr->ystep << "\t" << itr->zstep << "\t";
			std::cout << itr->xbigmax << "\t" << itr->ybigmax << "\t" << itr->zbigmax << "\t";
			std::cout << itr->xbigstep << "\t" << itr->ybigstep << "\t" << itr->zbigstep;

			std::cout << std::endl;

			itr++;
		}
	}
#endif

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
			std::cout << "GEANT4> RETURNING PDG id: " << *it << std::endl;
			it++;
		}
	}
	tracking->SetParticlesToKeep(&keep_ids);

	if(*g4_debug)
	{
		std::cout << "GEANT4> keep stable particles:  " << *g4_keep_stable << std::endl;
		std::cout << "GEANT4> keep neutral particles: " << *g4_neutral     << std::endl;
	}
	tracking->SetKeepStableParticles(*g4_keep_stable);
	tracking->SetKeepNeutrals(*g4_neutral);
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

	//Set up the stacking action to kill non-desired particles
	runManager->SetUserAction(stack);

	//Added everything now set up the run manager
	runManager->Initialize();

	MaximumParticleID = 0;

	std::cout << std::flush;
}

extern "C" void g4_add_collimator(char* name, char* material, double* length, double* aperture, double* rotation, double* x_offset, double* y_offset, bool* onesided)
{
	//NOTE: The closed orbit offset (e.g. at TCTs due to the crossing angle) should be subtracted in the sixtrack part.
	//rcx(j)  = (xv(1,j)-torbx(ie))/1d3
	//rcy(j)  = (xv(2,j)-torby(ie))/1d3
	//Therefore we do not need to take it into account here...

//  keep 48 value in sync with mNameLen in common_modules.f90
	std::string CollimatorName = CleanFortranString(name, 48);
	std::string CollimatorMaterialName = CleanFortranString(material, 4);
	std::cout << "GEANT4> Adding \"" << CollimatorName << "\" with material \"" << CollimatorMaterialName << "\" and rotation \"" << *rotation << "\" and offset x: \"" << *x_offset << "\" y: \"" << *y_offset << "\" and length \"";
	std::cout << *length << "\"" << std::endl;

	G4double length_in = *length *CLHEP::m;
	G4double aperture_in = *aperture *CLHEP::m;
	G4double rotation_in = *rotation *CLHEP::rad;
	G4double offset_in = *x_offset *CLHEP::m;

	geometry->AddCollimator(CollimatorName, length_in, aperture_in, rotation_in, offset_in, CollimatorMaterialName, *onesided);
}

extern "C" void g4_terminate()
{
	std::cout << "GEANT4> EXIT GEANT4" << std::endl;

	if(runManager)
	{
		std::cout << "GEANT4> Deleting geant4 run manager" << std::endl;
		delete runManager;
	}
}

//Set up new collimator
//runManager->ReinitializeGeometry();
extern "C" void g4_set_collimator(char* name)
{
//  keep 48 value in sync with mNameLen in common_modules.f90
	std::string CollimatorName = CleanFortranString(name, 48);
	geometry->SetCollimator(CollimatorName);

	G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

	runManager->ReinitializeGeometry();

	std::cout << std::flush;
}

extern "C" void g4_add_particle(double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, int16_t* nzz, int16_t* naa, int16_t* nqq, double* mass, double* sigmv, int32_t* partID, int32_t* parentID, double* weight, double* sx, double* sy, double* sz)
{
//WARNING: at this stage in SixTrack the units have been converted to GeV, m, and rad!
//The particle energy input is the TOTAL energy
//mass (i.e. nucm) is already in MeV!

	double x_in = (*x) * CLHEP::m;
	double y_in = (*y) * CLHEP::m;

//We want px and py, not the angle!
//p_in is the TOTAL momentum of the particle

	double e_in = (*e) * CLHEP::GeV;
	double p_in = sqrt((e_in*e_in) - (*mass * *mass));

// x' = p_x / p_in
// -> p_x = p_in * x'
//p_in will now be in MeV, xp, yp will be in rad -> units are good!
	double px_in = p_in * (*xp);
	double py_in = p_in * (*yp);

// p_z^2 = p_in^2 - p_x^2 - p_y^2
	double p_z = sqrt( (p_in*p_in) - (px_in*px_in) - (py_in*py_in) );

	CollimationParticle in_particle;
	in_particle.x = x_in;
	in_particle.y = y_in;

	in_particle.px = px_in;
	in_particle.py = py_in;
	in_particle.pz = p_z;
	in_particle.p = p_in;

	in_particle.e = e_in;
	in_particle.pdgid = *pdgid;
	in_particle.z = *nzz;
	in_particle.a = *naa;
	in_particle.q = *nqq;
	in_particle.m = *mass;

	in_particle.id = *partID;
	in_particle.parent_id = *parentID;
	in_particle.weight = *weight;

	in_particle.sx = *sx;
	in_particle.sy = *sy;
	in_particle.sz = *sz;

//  This should now be in seconds
	in_particle.t = *sigmv;

	input_particles.push_back(in_particle);

	if(*partID > MaximumParticleID)
	{
		MaximumParticleID = *partID;
		std::cout << "GEANT4> WARNING: Particle ID entered (" << *partID << ") is greater than the expected maximum currently set: " << MaximumParticleID << std::endl;
	}

}

extern "C" void g4_collimate()
{
	output_particles.clear();
	//Update the gun with this particle's details
	for(size_t n=0; n < input_particles.size(); n++)
	{
		part->SetParticleDetails(input_particles.at(n).x, input_particles.at(n).y, input_particles.at(n).px, input_particles.at(n).py, input_particles.at(n).pz, input_particles.at(n).e, input_particles.at(n).p, input_particles.at(n).t, input_particles.at(n).pdgid, input_particles.at(n).q, input_particles.at(n).m, input_particles.at(n).sx, input_particles.at(n).sy, input_particles.at(n).sz,input_particles.at(n).weight);

		//Tell the "tracking" about the parent particle ID for tracking secondaries
		tracking->SetParticleID(input_particles.at(n).id);
		tracking->SetParentID(input_particles.at(n).parent_id);
		tracking->SetMaximumParticleID(MaximumParticleID);
		//Run!
		runManager->BeamOn(1);
		MaximumParticleID = tracking->GetMaximumParticleID();
	}
	input_particles.clear();
	std::cout << std::flush;
}

/**
* Here we put the particles back into sixtrack and set any flags if needed
*/
extern "C" void g4_collimate_return(int* j, double* x, double* y, double* xp, double* yp, double* e, int32_t* pdgid, double* m, int16_t* z, int16_t* a, int16_t* q, double* sigmv, int32_t* partID, int32_t* parentID, double* weight, int *part_hit, int *part_abs, double *part_impact, double *part_indiv, double *part_linteract, double* sx, double* sy, double* sz)
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

//Here the units have been converted back to GeV and m (in the tracking action)

*x  = output_particles.at(*j).x;
*y  = output_particles.at(*j).y;
double px = output_particles.at(*j).px;
double py = output_particles.at(*j).py;

//Remember, sixtrack xp, yp are p_x / p_total
*xp = output_particles.at(*j).px / output_particles.at(*j).p;
*yp = output_particles.at(*j).py / output_particles.at(*j).p;
*e  = output_particles.at(*j).e;
*pdgid  = output_particles.at(*j).pdgid;
*z = output_particles.at(*j).z;
*a  = output_particles.at(*j).a;
*q  = output_particles.at(*j).q;

//nucm is in MeV on both sides
*m  = output_particles.at(*j).m;

//Spin
*sx = output_particles.at(*j).sx;
*sy = output_particles.at(*j).sy;
*sz = output_particles.at(*j).sz;

//time, must be converted for using with sigmv (done on the Sixtrack fortran side)
*sigmv  = output_particles.at(*j).t;


*partID = output_particles.at(*j).id;
*parentID = output_particles.at(*j).parent_id;
*weight = output_particles.at(*j).weight;
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

	std::transform(sstring.begin(), sstring.end(), sstring.begin(), ::toupper);
	//Fortran string happy place
	return sstring;
}

extern "C" void g4_get_particle_count(int* g4_npart)
{
	*g4_npart = output_particles.size();
}

extern "C" void g4_get_maximum_particle_id(int32_t* m_id)
{
	*m_id = tracking->GetMaximumParticleID();
}

//Sets the maximum particle ID (which could have been increased by FLUKA or other codes
//These particles could then be lost before geant4 encounters them (so the maximum must be tracked)!
extern "C" void g4_set_maximum_particle_id(int32_t* m_id)
{
	tracking->SetMaximumParticleID(*m_id);
}

extern "C" void g4_collimation_clear()
{
	input_particles.clear();
	output_particles.clear();
}

extern "C" void g4_keep_id(int* id)
{
	keep_ids.insert(*id);
}

extern "C" void g4_add_edep(char* name_in, int* type, double* xstep, double* ystep, double* zstep, double* xmax, double* ymax, double* zmax, double* xbigstep, double* ybigstep, double* zbigstep, double* xbigmax, double* ybigmax, double* zbigmax)
{
	std::string name(name_in);
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);

	//find name in struct vector
	std::vector<CollimationEnergyDeposition>::iterator itr = std::find_if(EnergyDepositionConfiguration.begin(), EnergyDepositionConfiguration.end(), find_edep_name(name));

	if(itr != EnergyDepositionConfiguration.end())
	{
		//update
		//Check for non-zero values and update
		if(*type != 0)
		{
			itr->type = *type;
		}
		if(*xstep != 0)
		{
			itr->xstep = *xstep;
		}
		if(*ystep != 0)
		{
			itr->ystep = *ystep;
		}
		if(*zstep != 0)
		{
			itr->zstep = *zstep;
		}

		if(*xbigstep != 0)
		{
			itr->xbigstep = *xbigstep;
		}
		if(*ybigstep != 0)
		{
			itr->ybigstep = *ybigstep;
		}
		if(*zbigstep != 0)
		{
			itr->zbigstep = *zbigstep;
		}

		if(*xmax != 0)
		{
			itr->xmax = *xmax;
		}
		if(*ymax != 0)
		{
			itr->ymax = *ymax;
		}
		if(*zmax != 0)
		{
			itr->zmax = *zmax;
		}

		if(*xbigmax != 0)
		{
			itr->xbigmax = *xbigmax;
		}
		if(*ybigmax != 0)
		{
			itr->ybigmax = *ybigmax;
		}
		if(*zbigmax != 0)
		{
			itr->zbigmax = *zbigmax;
		}
	}
	else
	{
		//if not found, create
		CollimationEnergyDeposition tmp;

		//Remember values are in mm
		double DefaultSmallStep = 0.05;
		double DefaultBigStep   = 1.0;

		double DefaultSmallMax = 10.0;
		double DefaultBigMax   = 40.0;

		//For z, assume default zmax is CollimatorLength : use -1 as a flag for this

		tmp.name = name;
		tmp.type = 0;
		tmp.xstep = DefaultSmallStep;
		tmp.ystep = DefaultSmallStep;
		tmp.zstep = DefaultSmallStep;

		tmp.xmax = DefaultSmallMax;
		tmp.ymax = DefaultSmallMax;
		tmp.zmax = -1;

		tmp.xbigstep = DefaultBigStep;
		tmp.ybigstep = DefaultBigStep;
		tmp.zbigstep = DefaultBigStep;

		tmp.xbigmax = DefaultBigMax;
		tmp.ybigmax = DefaultBigMax;
		tmp.zbigmax = -1;
		//Check for non-zero values and update
		if(*type != 0)
		{
			tmp.type = *type;
		}
		if(*xstep != 0)
		{
			tmp.xstep = *xstep;
		}
		if(*ystep != 0)
		{
			tmp.ystep = *ystep;
		}
		if(*zstep != 0)
		{
			tmp.zstep = *zstep;
		}

		if(*xbigstep != 0)
		{
			tmp.xbigstep = *xbigstep;
		}
		if(*ybigstep != 0)
		{
			tmp.ybigstep = *ybigstep;
		}
		if(*zbigstep != 0)
		{
			tmp.zbigstep = *zbigstep;
		}

		if(*xmax != 0)
		{
			tmp.xmax = *xmax;
		}
		if(*ymax != 0)
		{
			tmp.ymax = *ymax;
		}
		if(*zmax != 0)
		{
			tmp.zmax = *zmax;
		}

		if(*xbigmax != 0)
		{
			tmp.xbigmax = *xbigmax;
		}
		if(*ybigmax != 0)
		{
			tmp.ybigmax = *ybigmax;
		}
		if(*zbigmax != 0)
		{
			tmp.zbigmax = *zbigmax;
		}

		EnergyDepositionConfiguration.push_back(tmp);
	}
}

