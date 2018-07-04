#include "Collimation_root.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>

CollimationRootOutput* CollimationLossOutput;
CollimationDBRootOutput* CollimationDBOutput;
CollimationFLUKARootOutput* CollimationFLUKAOutput;

/**
* Collimation start up and tree creation
* Need collimation losses
* Aperture losses
*/
extern "C" void CollimationRootInit()
{
	CollimationLossOutput = new CollimationRootOutput();
}

extern "C" void CollimationDBRootInit()
{
    CollimationDBOutput = new CollimationDBRootOutput();
}

extern "C" void CollimationFLUKARootInit()
{
    CollimationFLUKAOutput = new CollimationFLUKARootOutput();
}

//Just pass values forward to c++
extern "C" void CollimatorLossRootWrite(int icoll_in, char* db_name_in, int db_name_len, int impact_in, int absorbed_in, double caverage_in, double csigma_in, double length_in)
{
	CollimationLossOutput->CollimationLossRootOutputWrite(icoll_in, db_name_in, db_name_len, impact_in, absorbed_in, caverage_in, csigma_in, length_in);
}

extern "C" void SurvivalRootWrite(int nturn_in, int npart_in)
{
    CollimationLossOutput->SurvivalRootOutputWrite(nturn_in, npart_in);
}

extern "C" void CollimatorDatabaseRootWrite(int j, char* db_name_in, int db_name_len, char* db_material_in, int db_material_len, double db_nsig_in, double db_length_in, double db_rotation_in, double db_offset_in)
{
    CollimationDBOutput->CollimatorDatabaseRootOutputWrite(j, db_name_in, db_name_len, db_material_in, db_material_len, db_nsig_in, db_length_in, db_rotation_in, db_offset_in);
}

extern "C" void root_FLUKA_EnergyDeposition(int id_in, int nucleons_in, double energy_in)
{
    CollimationFLUKAOutput->CollimatorFLUKARootOutputWrite(id_in, nucleons_in, energy_in);
}

CollimationRootOutput::CollimationRootOutput()
{
//Tree stuff
	CollimatorLossTree = new TTree("CollimatorLoss","CollimatorLossTree");
	CollimatorLossTree->Branch("icoll",&icoll,"icoll/I");
	CollimatorLossTree->Branch("impact",&impact,"impact/I");
	CollimatorLossTree->Branch("absorbed",&absorbed,"absorbed/I");
	CollimatorLossTree->Branch("caverage",&caverage,"caverage/D");
	CollimatorLossTree->Branch("csigma",&csigma,"csigma/D");
	CollimatorLossTree->Branch("length",&length,"length/D");
	CollimatorLossTree->Branch("db_name",name,"db_name[49]/C");

	CollimationSurvivalTree = new TTree("CollimationSurvival","CollimationSurvivalTree");
	CollimationSurvivalTree->Branch("nturn",&nturn,"nturn/I");
	CollimationSurvivalTree->Branch("npart",&npart,"npart/I");
}

void CollimationRootOutput::CollimationLossRootOutputWrite(int icoll_in, char* db_name_in, int db_name_len, int impact_in, int absorbed_in, double caverage_in, double csigma_in, double length_in)
{
	icoll = icoll_in;
	impact = impact_in;
	absorbed = absorbed_in;
	caverage = caverage_in;
	csigma = csigma_in;
	length = length_in;

	std::string temp_name(db_name_in);
	std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::toupper);
	strncpy(name,temp_name.substr(0,db_name_len).c_str(),db_name_len);

	CollimatorLossTree->Fill();
}


void CollimationRootOutput::SurvivalRootOutputWrite(int nturn_in, int npart_in)
{
	nturn = nturn_in;
	npart = npart_in;
	CollimationSurvivalTree->Fill();
}

/**
* Collimator database functions
*/

CollimationDBRootOutput::CollimationDBRootOutput()
{
	CollimatorDatabaseTree = new TTree("CollimatorDatabase","CollimatorDatabaseTree");
	CollimatorDatabaseTree->Branch("icoll",&icoll,"icoll/I");
	CollimatorDatabaseTree->Branch("db_nsig",&db_nsig,"db_nsig/D");
	CollimatorDatabaseTree->Branch("db_length",&db_length,"db_length/D");
	CollimatorDatabaseTree->Branch("db_rotation",&db_rotation,"db_rotation/D");
	CollimatorDatabaseTree->Branch("db_offset",&db_offset,"db_offset/D");
	CollimatorDatabaseTree->Branch("db_name",db_name,"db_name[49]/C");
	CollimatorDatabaseTree->Branch("db_material",db_material,"db_material[5]/C");
}

void CollimationDBRootOutput::CollimatorDatabaseRootOutputWrite(int j, char* db_name_in, int db_name_len, char* db_material_in, int db_material_len, double db_nsig_in, double db_length_in, double db_rotation_in, double db_offset_in)
{
	icoll = j;

	std::string temp_name(db_name_in);
	std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::toupper);
	strncpy(db_name,temp_name.substr(0,db_name_len).c_str(),db_name_len);

	std::string temp_material(db_material_in);
	strncpy(db_material,temp_material.substr(0,db_material_len).c_str(),db_material_len);

	db_nsig = db_nsig_in;
	db_length = db_length_in;
	db_rotation = db_rotation_in;
	db_offset = db_offset_in;
//	std::cout << "\"" << db_name << "\"\t\"" << db_material << "\"\tis at " << db_nsig << " sigma" << std::endl;
//	std::cout << db_name_len << "\t" << db_material_len << std::endl;
	CollimatorDatabaseTree->Fill();
}

/**
* Collimator FLUKA output 
*/

CollimationFLUKARootOutput::CollimationFLUKARootOutput()
{
	CollimationFLUKATree = new TTree("CollimatorFLUKA","CollimationFLUKATree");
	CollimationFLUKATree->Branch("id",&id,"id/I");
	CollimationFLUKATree->Branch("nucleons",&nucleons,"nucleons/I");
	CollimationFLUKATree->Branch("energy",&energy,"energy/D");
}

void CollimationFLUKARootOutput::CollimatorFLUKARootOutputWrite(int id_in, int nucleons_in, double energy_in)
{
	id = id_in;

	nucleons = nucleons_in;
	energy = energy_in;

	CollimationFLUKATree->Fill();
}

