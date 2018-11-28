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

#include "CollimationGeometry.h"

//This has the G4 collimator geometry construction
//Size, material, length, gap sizes, tilts.
/*
void BuildCollimator(double material, double length, double gap)
{

}
*/
G4VPhysicalVolume* CollimationGeometry::Construct()
{
	if(Assembled == true)
	{

		Assembled = false;
	}
	//Change this to true if there are future geometry issues
	bool DoDebug = false;
	bool OverlapCheck = false;
	if(DoDebug)
	{
		OverlapCheck = true;
		std::cout << "GAP :" << ThisCollimatorJawHalfGap << std::endl;
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
	world_log = new G4LogicalVolume(world_box, Vacuum, "world_log");
	Jaw1_log  = new G4LogicalVolume(Jaw1,ThisCollimatorJawMaterial,"Jaw1_log");
	Jaw2_log  = new G4LogicalVolume(Jaw2,ThisCollimatorJawMaterial,"Jaw2_log");

	world_phys = new G4PVPlacement(0,G4ThreeVector(), world_log, "world", 0, false, 0);
	jaw1_phys = new G4PVPlacement(0, G4ThreeVector( (jaw_x)+HalfGap, 0, 0.5*Length), Jaw1_log, "jaw1", world_log, false, 0, OverlapCheck);
	jaw2_phys = new G4PVPlacement(0, G4ThreeVector(-(jaw_x)-HalfGap, 0, 0.5*Length), Jaw2_log, "jaw2", world_log, false, 0, OverlapCheck);

	if(DoDebug)
	{
		G4cout << "Adding new collimator with name: " << ThisCollimatorName << " and material " << ThisCollimatorJawMaterial->GetName() << std::endl;
		G4cout << "Total Jaw Length: " << Length/CLHEP::m <<  "m" << G4endl;
		G4cout << "Jaw Rotation: " << ThisCollimatorJawRotation/CLHEP::rad <<  "rad" << G4endl;
		G4cout << "+ve Jaw position: " << ((jaw_x)+HalfGap)/CLHEP::m << "m" << G4endl;
		G4cout << "-ve Jaw position: " << (-(jaw_x)-HalfGap)/CLHEP::m << "m" << G4endl;
		G4cout << "Jaw Gap: " << (HalfGap)/CLHEP::mm << "mm" << G4endl;
//		G4cout << "t1: " << position1 << G4endl; 
//		G4cout << "t1: " << rotation << G4endl; 
//		G4cout << "t2: " << position2 << G4endl; 
	}

	Assembled = true;
	return world_phys;
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
		std::cerr << "ERROR: Multiple definitions of collimator: \"" << name << "\"" << std::endl;
		exit(EXIT_FAILURE);
	}
}

void CollimationGeometry::SetCollimator(std::string CollimatorName)
{
//	std::cout << "Looking for collimator \"" << CollimatorName << "\"" << std::endl;
	std::map<std::string,CollimatorSettings>::iterator CollimatorKey_itr = CollimatorKeyMap.find(CollimatorName);
	if(CollimatorKey_itr == CollimatorKeyMap.end())
	{
		std::cerr << "ERROR: Collimator \"" <<  CollimatorName << "\" not found!" << std::endl;
		exit(EXIT_FAILURE);
	}
	ThisCollimatorJawLength = CollimatorKey_itr->second.CollimatorJawLength;
	ThisCollimatorJawHalfGap = CollimatorKey_itr->second.CollimatorJawHalfGap;
	ThisCollimatorJawRotation = CollimatorKey_itr->second.CollimatorJawRotation;
	ThisCollimatorJawOffset = CollimatorKey_itr->second.CollimatorJawOffset;
	ThisCollimatorJawMaterial = CollimatorKey_itr->second.CollimatorJawMaterial;
	ThisCollimatorName = CollimatorKey_itr->second.CollimatorName;
}
