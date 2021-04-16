#include <algorithm>
#include <string>

#include "CollimationMaterials.h"

/*
BE
AL
CU
W
PB
C
C2
MoGR
CuCD
Mo
Glid
Iner
*/

CollimationMaterials::CollimationMaterials()
{
	G4NistManager* NManager = G4NistManager::Instance();

	//Elements needed for composite materials
	G4Material* B = NManager->FindOrBuildMaterial("G4_B");
	G4Material* Be = NManager->FindOrBuildMaterial("G4_Be");
	G4Material* C = NManager->FindOrBuildMaterial("G4_C");
	G4Material* O = NManager->FindOrBuildMaterial("G4_O");
	G4Material* Al = NManager->FindOrBuildMaterial("G4_Al");
	G4Material* Cu = NManager->FindOrBuildMaterial("G4_Cu");
	G4Material* Ni = NManager->FindOrBuildMaterial("G4_Ni");
	G4Material* Mo = NManager->FindOrBuildMaterial("G4_Mo");
	G4Material* W = NManager->FindOrBuildMaterial("G4_W");
	G4Material* Pb = NManager->FindOrBuildMaterial("G4_Pb");
	G4Material* Si = NManager->FindOrBuildMaterial("G4_Si");


	//Only needed for pure elemental collimators
	AddMaterial("BE", Be);
//	AddMaterial("C", C);
	AddMaterial("O", O);
	AddMaterial("AL", Al);
	AddMaterial("CU", Cu);
	AddMaterial("NI", Ni);
	AddMaterial("MO", Mo);
	AddMaterial("W", W);
	AddMaterial("PB", Pb);
	AddMaterial("SI", Si);

// Remember:
//  void AddMaterial(G4Material* material,                        //the material
//                   G4double   fraction);                        //fractionOfMass
// Mass fractions!
//
	G4Material* AC150K = new G4Material("AC150K", 1.670*CLHEP::g/CLHEP::cm3,1);
	AC150K->AddMaterial(C,1.0);
	AddMaterial("C", AC150K);

	//Mo Graphite - fractions from https://twiki.cern.ch/twiki/pub/LHCAtHome/SixTrackCollimatVer/material_test_2015-03-30_corretto.xlsx
	G4Material* MoGr = new G4Material("MoGr", 2.5*CLHEP::g/CLHEP::cm3,2);
	MoGr->AddMaterial(Mo,0.137);
	MoGr->AddMaterial(C,0.863);
	AddMaterial("MOGR", MoGr);

	//Glidcop
	G4Material* Glid = new G4Material("Glid", 8.930*CLHEP::g/CLHEP::cm3,3);
	G4double Al_M = 0.4 * Al->GetA();
	G4double O_M = 0.6 * O->GetA();
	Glid->AddMaterial(Cu,0.9972);
	Glid->AddMaterial(Al,0.0028 * Al_M / (Al_M + O_M));
	Glid->AddMaterial(O,0.0028 * O_M / (Al_M + O_M));
	AddMaterial("GLID", Glid);

	//INERMET 180
	G4Material* Iner = new G4Material("Iner", 18.060*CLHEP::g/CLHEP::cm3,3);
	Iner->AddMaterial(W,0.95);
	Iner->AddMaterial(Ni,0.035);
	Iner->AddMaterial(Cu,0.015);
	AddMaterial("INER", Iner);

	//Copper diamond - fractions from https://twiki.cern.ch/twiki/pub/LHCAtHome/SixTrackCollimatVer/material_test_2015-03-30_corretto.xlsx
	//https://cds.cern.ch/record/2112203/files/tho4ab03.pdf
	G4Material* CuCD = new G4Material("CuCD", 5.4*CLHEP::g/CLHEP::cm3,3);
	CuCD->AddMaterial(Cu,0.647);
	CuCD->AddMaterial(C,0.349);
	CuCD->AddMaterial(B,0.004);
	AddMaterial("CUCD", CuCD);
}


void CollimationMaterials::AddMaterial(std::string MaterialName, G4Material* Material)
{
	std::pair<std::string, G4Material*> MaterialMapKey;
	MaterialMapKey = std::make_pair(MaterialName, Material);
	std::pair<std::map<std::string,G4Material*>::iterator,bool> check = CollimationMaterialMap.insert(MaterialMapKey);
	if(!check.second)
	{
		std::cerr << "GEANT4> ERROR: Failed to insert entry into G4Collimation material storage!: " << MaterialName << std::endl;
		exit(EXIT_FAILURE);

	}
}

G4Material* CollimationMaterials::GetMaterial(std::string MaterialName)
{
	std::transform(MaterialName.begin(), MaterialName.end(), MaterialName.begin(), ::toupper);
	std::map<std::string, G4Material*>::const_iterator CollimationMaterialMap_itr;
	CollimationMaterialMap_itr = CollimationMaterialMap.find(MaterialName);
	if(CollimationMaterialMap_itr!=CollimationMaterialMap.end())
	{
		return CollimationMaterialMap_itr->second;
	}
	else
	{
		G4cout << "ERROR: Material \"" << MaterialName << "\" was requested but was not found!" << G4endl;
		exit(EXIT_FAILURE);
	}
}
