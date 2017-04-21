#ifndef _CollimationGeometry_h
#define _CollimationGeometry_h

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include <string>

class CollimationGeometry : public G4VUserDetectorConstruction
{
	public:

	CollimationGeometry(std::string name, double length, double gap, double rotation, double offset, G4Material* Material)
    : CollimatorJawLength(length), CollimatorJawHalfGap(gap*0.5), CollimatorJawRotation(rotation), CollimatorJawOffset(offset), CollimatorJawMaterial(Material), CollimatorName(name) {};

	G4VPhysicalVolume* Construct();

	G4double GetLength();
	G4double GetHalfGap();

	private:

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


	G4Box* world_box;
	G4Box* Jaw1;
	G4Box* Jaw2;
	G4LogicalVolume* world_log;
	G4LogicalVolume* Jaw1_log;
	G4LogicalVolume* Jaw2_log;
	G4PVPlacement* world_phys;
	G4VPhysicalVolume* jaw1_phys;
	G4VPhysicalVolume* jaw2_phys;

};

#endif

