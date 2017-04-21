#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

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
	//Change this to true if there are future geometry issues
	bool DoDebug = false;
	bool OverlapCheck = false;
	if(DoDebug)
	{
		OverlapCheck = true;
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
	G4double HalfGap = CollimatorJawHalfGap;
	G4double Length = CollimatorJawLength;

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
	Jaw1_log  = new G4LogicalVolume(Jaw1,CollimatorJawMaterial,"Jaw1_log");
	Jaw2_log  = new G4LogicalVolume(Jaw2,CollimatorJawMaterial,"Jaw2_log");

	world_phys = new G4PVPlacement(0,G4ThreeVector(), world_log, "world", 0, false, 0);
/*
	G4RotationMatrix rotation  = G4RotationMatrix();
	rotation.rotateZ(CollimatorJawRotation);
	G4double JawPosition = jaw_x + HalfGap;
	G4double pi = 4.0*atan2(1.0,1.0);

	//Placement point: rotate around x,y and place a distance of the jaw gap away. Then rotate the jaw face!
	G4ThreeVector position1 = G4ThreeVector(JawPosition*cos(CollimatorJawRotation), JawPosition*sin(CollimatorJawRotation), 0.5 * Length);
	G4ThreeVector position2 = G4ThreeVector(JawPosition*cos(CollimatorJawRotation+pi), JawPosition*sin(CollimatorJawRotation+pi), 0.5 * Length);
	const G4Transform3D transform1 = G4Transform3D(rotation, position1);
	const G4Transform3D transform2 = G4Transform3D(rotation, position2);

	jaw1_phys = new G4PVPlacement(transform1, Jaw1_log, "jaw1", world_log, false, 0, OverlapCheck);
	jaw2_phys = new G4PVPlacement(transform2, Jaw2_log, "jaw2", world_log, false, 0, OverlapCheck);
*/
	jaw1_phys = new G4PVPlacement(0, G4ThreeVector( (jaw_x)+HalfGap, 0, 0.5*Length), Jaw1_log, "jaw1", world_log, false, 0, OverlapCheck);
	jaw2_phys = new G4PVPlacement(0, G4ThreeVector(-(jaw_x)-HalfGap, 0, 0.5*Length), Jaw2_log, "jaw2", world_log, false, 0, OverlapCheck);

	if(DoDebug)
	{
		G4cout << "Adding new collimator with name: " << CollimatorName << " and material " << CollimatorJawMaterial->GetName() << std::endl;
		G4cout << "Total Jaw Length: " << Length/CLHEP::m <<  "m" << G4endl; 
		G4cout << "Jaw Rotation: " << CollimatorJawRotation/CLHEP::rad <<  "rad" << G4endl; 
		G4cout << "+ve Jaw position: " << ((jaw_x)+HalfGap)/CLHEP::m << "m" << G4endl; 
		G4cout << "-ve Jaw position: " << (-(jaw_x)-HalfGap)/CLHEP::m << "m" << G4endl; 
		G4cout << "Jaw Gap: " << (HalfGap)/CLHEP::mm << "mm" << G4endl; 
//		G4cout << "t1: " << position1 << G4endl; 
//		G4cout << "t1: " << rotation << G4endl; 
//		G4cout << "t2: " << position2 << G4endl; 
	}
	return world_phys;
}

G4double CollimationGeometry::GetLength()
{
	return CollimatorJawLength;
}

G4double CollimationGeometry::GetHalfGap()
{
	return CollimatorJawHalfGap;
}
