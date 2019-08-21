#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "CollimationGeometry.h"
#include "CollimationJawSD.h"
#include "CollimationEventAction.h"

//This has the G4 collimator geometry construction
//Size, material, length, gap sizes, tilts.
G4VPhysicalVolume* CollimationGeometry::Construct()
{
	if(Assembled == true)
	{
		Assembled = false;
	}
	//Change this to true if there are future geometry issues
	bool OverlapCheck = false;
	if(DoDebug)
	{
		OverlapCheck = true;
		std::cout << "GEANT4> ThisCollimatorJawHalfGap :" << ThisCollimatorJawHalfGap << std::endl;
	}

	//Make materials
	G4NistManager* NManager = G4NistManager::Instance();
	G4Material* Vacuum = NManager->FindOrBuildMaterial("G4_Galactic");

	/*
	* Configure the master "world" box
	* "Note that the G4Box constructor takes as arguments the halfs of the total box size."
	*/
	//convert input into G4 units
	const G4double jaw_x = 50*CLHEP::m / 2.0;
	const G4double jaw_y = 50*CLHEP::m / 2.0;

	//The number from SixTrack is the full gap size
	G4double HalfGap = ThisCollimatorJawHalfGap;
	G4double Length = ThisCollimatorJawLength;

	//We are going to place the gun at zero to keep things simple
	//The box spans from -Jaw length to +Jaw length
	G4double world_x = 100*CLHEP::m;
	G4double world_y = 100*CLHEP::m;
	G4double world_z = Length;

	//Make the world box
	world_box = new G4Box("world_box", world_x, world_y, world_z);

	//Make the left/right jaws
	Jaw1 = new G4Box("Jaw1_box", jaw_x, jaw_y, 0.5*Length);
	Jaw2 = new G4Box("Jaw2_box", jaw_x, jaw_y, 0.5*Length);

	//Make the logical volumes, and attach materials
	world_log = new G4LogicalVolume(world_box, Vacuum,                    "world_log");
	Jaw1_log  = new G4LogicalVolume(Jaw1,      ThisCollimatorJawMaterial, "Jaw1_log");
	Jaw2_log  = new G4LogicalVolume(Jaw2,      ThisCollimatorJawMaterial, "Jaw2_log");

	world_phys = new G4PVPlacement(0, G4ThreeVector(),                                world_log, "world", 0,         false, 0);
	jaw1_phys  = new G4PVPlacement(0, G4ThreeVector( (jaw_x)+HalfGap, 0, 0.5*Length), Jaw1_log,  "jaw1",  world_log, false, 0, OverlapCheck);
	jaw2_phys  = new G4PVPlacement(0, G4ThreeVector(-(jaw_x)-HalfGap, 0, 0.5*Length), Jaw2_log,  "jaw2",  world_log, false, 0, OverlapCheck);

	if(DoDebug)
	{
		G4int OldPrecision = G4cout.precision(9);
		G4cout << "GEANT4> Adding new collimator with name: " << ThisCollimatorName << " and material " << ThisCollimatorJawMaterial->GetName() << std::endl;
		G4cout << "GEANT4> Total Jaw Length: " << Length/CLHEP::m <<  "m" << G4endl;
		G4cout << "GEANT4> Jaw Rotation: " << ThisCollimatorJawRotation/CLHEP::rad <<  "rad" << G4endl;
		G4cout << "GEANT4> +ve Jaw position: " << ((jaw_x)+HalfGap)/CLHEP::m << "m" << G4endl;
		G4cout << "GEANT4> -ve Jaw position: " << (-(jaw_x)-HalfGap)/CLHEP::m << "m" << G4endl;
		G4cout << "GEANT4> Jaw Half Gap: " << (HalfGap)/CLHEP::mm << "mm" << G4endl;
		G4cout.precision(OldPrecision);
	}

	Assembled = true;
	return world_phys;
}


CollimationGeometry::~CollimationGeometry()
{

	if(world_phys)
		delete world_phys;
	if(jaw1_phys)
		delete jaw1_phys;
	if(jaw2_phys)
		delete jaw2_phys;


	if(world_log)
		delete world_log;
	if(Jaw1_log)
		delete Jaw1_log;
	if(Jaw2_log)
		delete Jaw2_log;

	if(Jaw1)
		delete Jaw1;
	if(Jaw2)
		delete Jaw2;

	if(world_box)
		delete world_box;
}

void CollimationGeometry::ConstructSDandField()
{
#ifdef USE_ROOT_FLAG
	if(DoEnergyDeposition)
	{
	G4SDManager* sdm = G4SDManager::GetSDMpointer();

	std::string Left  = ThisCollimatorName+"_LeftJaw";
	std::string Right = ThisCollimatorName+"_RightJaw";

	//We need to check if this class already exists, and if so, re-use it
	G4VSensitiveDetector* LeftJaw_SD_t = sdm->FindSensitiveDetector(Left, false);
	if(!LeftJaw_SD_t)
	{
		LeftJaw_SD  = new CollimationJawSD(Left, Left);
	}
	else
	{
		LeftJaw_SD = dynamic_cast<CollimationJawSD*>(LeftJaw_SD_t);
	}

	G4VSensitiveDetector* RightJaw_SD_t = sdm->FindSensitiveDetector(Right, false);
	if(!RightJaw_SD_t)
	{
		RightJaw_SD = new CollimationJawSD(Right, Right);
	}
	else
	{
		RightJaw_SD = dynamic_cast<CollimationJawSD*>(RightJaw_SD_t);
	}

	sdm->AddNewDetector(LeftJaw_SD);
	sdm->AddNewDetector(RightJaw_SD);

//	SetSensitiveDetector(G4LogicalVolume* logVol, G4VSensitiveDetector* aSD)

	SetSensitiveDetector("Jaw1_log", LeftJaw_SD);
	SetSensitiveDetector("Jaw2_log", RightJaw_SD);
	}
#endif
}


void CollimationGeometry::AddCollimator(std::string name, double length, double gap, double rotation, double offset, std::string Material)
{
	CollimatorSettings NewCollimator;
	NewCollimator.CollimatorJawLength = length;
	NewCollimator.CollimatorJawHalfGap = gap * 0.5;
	NewCollimator.CollimatorJawRotation = rotation;
	NewCollimator.CollimatorJawOffset = offset;
	G4Material* JawMaterial = Mmap->GetMaterial(Material);
	NewCollimator.CollimatorJawMaterial = JawMaterial;
	NewCollimator.CollimatorName = name;

	std::pair<std::string,CollimatorSettings> CollimatorKey;
	CollimatorKey = std::make_pair(name,NewCollimator);

	std::pair<std::map< std::string,CollimatorSettings >::const_iterator,bool > check;
	check = CollimatorKeyMap.insert(CollimatorKey);

	if(check.second == false)
	{
		std::cerr << "GEANT4> ERROR: Multiple definitions of collimator: \"" << name << "\"" << std::endl;
		exit(EXIT_FAILURE);
	}
}

void CollimationGeometry::SetCollimator(std::string CollimatorName)
{
//	std::cout << "Looking for collimator \"" << CollimatorName << "\"" << std::endl;
	std::map<std::string,CollimatorSettings>::iterator CollimatorKey_itr = CollimatorKeyMap.find(CollimatorName);
	if(CollimatorKey_itr == CollimatorKeyMap.end())
	{
		std::cerr << "GEANT4> ERROR: Collimator \"" <<  CollimatorName << "\" not found!" << std::endl;
		exit(EXIT_FAILURE);
	}
	ThisCollimatorJawLength = CollimatorKey_itr->second.CollimatorJawLength;
	ThisCollimatorJawHalfGap = CollimatorKey_itr->second.CollimatorJawHalfGap;
	ThisCollimatorJawRotation = CollimatorKey_itr->second.CollimatorJawRotation;
	ThisCollimatorJawOffset = CollimatorKey_itr->second.CollimatorJawOffset;
	ThisCollimatorJawMaterial = CollimatorKey_itr->second.CollimatorJawMaterial;
	ThisCollimatorName = CollimatorKey_itr->second.CollimatorName;

	//Pass this info to the event action
	((CollimationEventAction*)G4RunManager::GetRunManager()->GetUserEventAction())->SetCollimator(CollimatorName, ThisCollimatorJawLength, ThisCollimatorJawHalfGap);
}

void CollimationGeometry::SetDebug(bool flag)
{
	DoDebug = flag;
}

double CollimationGeometry::GetCurrentCollimatorLength()
{
	return ThisCollimatorJawLength;
}
