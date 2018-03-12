#include "Collimation_root.h"

#include "TTree.h"
#include "TROOT.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>

Int_t icoll;
Int_t impact;
Int_t absorbed;
Double_t caverage;
Double_t csigma;
Double_t length;

Int_t nturn;
Int_t npart;

//Collimator database variables
Char_t db_name[49];
Char_t db_material[5];
Double_t db_nsig;
Double_t db_length;
Double_t db_rotation;
Double_t db_offset;


TTree *CollimatorLossTree;
TTree *CollimationSurvivalTree;
TTree *CollimatorDatabaseTree;

/**
* Collimation start up and tree creation
* Need collimation losses
* Aperture losses
*/
extern "C" void CollimationRootInit()
{
//Tree stuff
    CollimatorLossTree = new TTree("CollimatorLoss","CollimatorLossTree");
    CollimatorLossTree->Branch("icoll",&icoll,"icoll/I");
    CollimatorLossTree->Branch("impact",&impact,"impact/I");
    CollimatorLossTree->Branch("absorbed",&absorbed,"absorbed/I");
    CollimatorLossTree->Branch("caverage",&caverage,"caverage/D");
    CollimatorLossTree->Branch("csigma",&csigma,"csigma/D");
    CollimatorLossTree->Branch("length",&length,"length/D");
    CollimatorLossTree->Branch("db_name",db_name,"db_name[49]/C");

    CollimationSurvivalTree = new TTree("CollimationSurvival","CollimationSurvivalTree");
    CollimationSurvivalTree->Branch("nturn",&nturn,"nturn/I");
    CollimationSurvivalTree->Branch("npart",&npart,"npart/I");

}

extern "C" void CollimationDBRootInit()
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

extern "C" void CollimatorLossRootWrite(int icoll_in, char* db_name_in, int db_name_len, int impact_in, int absorbed_in, double caverage_in, double csigma_in, double length_in)
{
    icoll = icoll_in;
    impact = impact_in;
    absorbed = absorbed_in;
    caverage = caverage_in;
    csigma = csigma_in;
    length = length_in;

    std::string temp_name(db_name_in);
    std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::toupper);
    strncpy(db_name,temp_name.substr(0,db_name_len).c_str(),db_name_len);

    CollimatorLossTree->Fill();
}

extern "C" void SurvivalRootWrite(int nturn_in, int npart_in)
{
    nturn = nturn_in;
    npart = npart_in;
    CollimationSurvivalTree->Fill();
}

extern "C" void CollimatorDatabaseRootWrite(int j, char* db_name_in, int db_name_len, char* db_material_in, int db_material_len, double db_nsig_in, double db_length_in, double db_rotation_in, double db_offset_in)
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
//    std::cout << "\"" << db_name << "\"\t\"" << db_material << "\"\tis at " << db_nsig << " sigma" << std::endl;
//    std::cout << db_name_len << "\t" << db_material_len << std::endl;
    CollimatorDatabaseTree->Fill();
}

