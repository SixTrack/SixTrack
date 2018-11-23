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


	//Only needed for pure elemental collimators
	AddMaterial("Be", Be);
	AddMaterial("C", C);
	AddMaterial("O", O);
	AddMaterial("Al", Al);
	AddMaterial("Cu", Cu);
	AddMaterial("Ni", Ni);
	AddMaterial("Mo", Mo);
	AddMaterial("W", W);
	AddMaterial("Pb", Pb);

	AddMaterial("CU", Cu);

	G4Material* AC150K = new G4Material("AC150K", 1.650*CLHEP::g/CLHEP::cm3,1);
	AC150K->AddMaterial(C,1.0);
	AddMaterial("C", AC150K);

	//Mo Graphite - FIXME fractions
	G4Material* MoGr = new G4Material("MoGr", 2.5*CLHEP::g/CLHEP::cm3,2);
	MoGr->AddMaterial(Mo,0.5);
	MoGr->AddMaterial(C,0.5);
	AddMaterial("MoGr", MoGr);
	AddMaterial("MoGR", MoGr);

	//Glidcop
	G4Material* Glid = new G4Material("Glid", 8.930*CLHEP::g/CLHEP::cm3,3);
	G4double Al_M = 0.4 * Al->GetA();
	G4double O_M = 0.6 * O->GetA();
	Glid->AddMaterial(Cu,0.9972);
	Glid->AddMaterial(Al,0.0028 * Al_M / (Al_M + O_M));
	Glid->AddMaterial(O,0.0028 * O_M / (Al_M + O_M));
	AddMaterial("Glid", Glid);

	//INERMET 180
	G4Material* Iner = new G4Material("Iner", 18.060*CLHEP::g/CLHEP::cm3,3);
	Iner->AddMaterial(W,0.95);
	Iner->AddMaterial(Ni,0.035);
	Iner->AddMaterial(Cu,0.015);
	AddMaterial("Iner", Iner);

	//Copper diamond - FIXME fractions
	//https://cds.cern.ch/record/2112203/files/tho4ab03.pdf
	G4Material* CuCD = new G4Material("CuCD", 5.4*CLHEP::g/CLHEP::cm3,3);
	CuCD->AddMaterial(Cu,0.39);
	CuCD->AddMaterial(C,0.6);
	CuCD->AddMaterial(B,0.01);
	AddMaterial("CuCD", CuCD);
}


void CollimationMaterials::AddMaterial(std::string MaterialName, G4Material* Material)
{
	std::pair<std::string, G4Material*> MaterialMapKey;
	MaterialMapKey = std::make_pair(MaterialName, Material);
	CollimationMaterialMap.insert(MaterialMapKey);
}

G4Material* CollimationMaterials::GetMaterial(std::string MaterialName)
{
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
