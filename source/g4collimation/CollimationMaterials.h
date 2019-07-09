#ifndef _CollimationMaterials_h_
#define _CollimationMaterials_h_ 1

#include <map>
#include <string>
#include "G4Material.hh"
#include "G4NistManager.hh"
class CollimationMaterials
{
	public:

	CollimationMaterials();

	//Add a created material
	void AddMaterial(std::string MaterialName, G4Material* Material);

	//Get a material
	G4Material* GetMaterial(std::string);

	private:
	std::map<std::string,G4Material*> CollimationMaterialMap;
};

#endif

