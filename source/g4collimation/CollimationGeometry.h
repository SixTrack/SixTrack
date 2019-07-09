#ifndef _CollimationGeometry_h
#define _CollimationGeometry_h

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "CollimationMaterials.h"

#include <string>
#include <map>

class CollimationGeometry : public G4VUserDetectorConstruction
{
	public:

	CollimationGeometry()
	{
		Mmap = new CollimationMaterials();

		Assembled = false;
		DoDebug = false;
		world_box = nullptr;
		Jaw1 = nullptr;
		Jaw2 = nullptr;
		world_log = nullptr;
		Jaw1_log = nullptr;
		Jaw2_log = nullptr;
		world_phys = nullptr;
		jaw1_phys = nullptr;
		jaw2_phys = nullptr;

		ThisCollimatorJawLength = 1.0;
		ThisCollimatorJawHalfGap = 1.0;
		ThisCollimatorJawRotation = 0.0;
		ThisCollimatorJawOffset = 0.0;
		ThisCollimatorJawMaterial = Mmap->GetMaterial("Iner");
		ThisCollimatorName = "dummy";
	};

	G4VPhysicalVolume* Construct();
	void SetCollimator(std::string);
	void AddCollimator(std::string name, double length, double gap, double rotation, double offset, std::string Material);
	void SetDebug(bool flag);


	private:

	struct CollimatorSettings
	{
		//The total length of the jaw in G4 units
		G4double CollimatorJawLength;

		//Half the total jaw gap in G4 units
		G4double CollimatorJawHalfGap;

		//Collimator jaw rotation
		G4double CollimatorJawRotation;

		//Collimator jaw offset
		G4double CollimatorJawOffset;

		//A pointer to a Geant4 material class containing the collimator jaw material
		G4Material* CollimatorJawMaterial;

		//The name of this collimator
		std::string CollimatorName;
	};

	std::map<std::string,CollimatorSettings> CollimatorKeyMap;

	bool Assembled;
	bool DoDebug;

	G4Box* world_box;
	G4Box* Jaw1;
	G4Box* Jaw2;
	G4LogicalVolume* world_log;
	G4LogicalVolume* Jaw1_log;
	G4LogicalVolume* Jaw2_log;
	G4PVPlacement* world_phys;
	G4VPhysicalVolume* jaw1_phys;
	G4VPhysicalVolume* jaw2_phys;

	CollimationMaterials* Mmap;

	G4double ThisCollimatorJawLength;
	G4double ThisCollimatorJawHalfGap;
	G4double ThisCollimatorJawRotation;
	G4double ThisCollimatorJawOffset;
	G4Material* ThisCollimatorJawMaterial;
	std::string ThisCollimatorName;
};

#endif

