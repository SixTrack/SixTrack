#ifndef Collimation_root_h
#define Collimation_root_h

#include "TTree.h"
#include "TROOT.h"

extern "C" void CollimationRootInit();
extern "C" void CollimationDBRootInit();
extern "C" void CollimationFLUKARootInit();
extern "C" void CollimationEnergyRootInit();

extern "C" void CollimatorLossRootWrite(int, char*, int, int, int, double, double, double);
extern "C" void SurvivalRootWrite(int, int);
extern "C" void CollimatorDatabaseRootWrite(int, char*, int, char*, int, double, double, double, double);
extern "C" void root_EnergyDeposition(int, int16_t, double);
extern "C" void root_FLUKA_Names(int, char*, int, int);


/**
* This class outputs the particle losses onto collimators
* It also prints the number of surviving particles each turn
*/
class CollimationRootOutput
{
public:

    CollimationRootOutput();
    void CollimationLossRootOutputWrite(int icoll_in, char* db_name_in, int db_name_len, int impact_in, int absorbed_in, double caverage_in, double csigma_in, double length_in);
    void SurvivalRootOutputWrite(int nturn_in, int npart_in);

private:

TTree *CollimatorLossTree;
TTree *CollimationSurvivalTree;

Char_t name[49];
Int_t icoll;
Int_t impact;
Int_t absorbed;
Double_t caverage;
Double_t csigma;
Double_t length;

Int_t nturn;
Int_t npart;

};

/**
* This class outputs the collimator database settings to root
*
*/
class CollimationDBRootOutput
{
public:

    CollimationDBRootOutput();
    void CollimatorDatabaseRootOutputWrite(int j, char* db_name_in, int db_name_len, char* db_material_in, int db_material_len, double db_nsig_in, double db_length_in, double db_rotation_in, double db_offset_in);

private:
TTree *CollimatorDatabaseTree;
//Collimator database variables
Char_t db_name[49];
Char_t db_material[5];
Double_t db_nsig;
Double_t db_length;
Double_t db_rotation;
Double_t db_offset;
Int_t icoll;

};

/**
* This class outputs the energy and nucleons dumped in collimators due to FLUKA/Geant4/etc
*/
class CollimationEnergyRootOutput
{
public:

    CollimationEnergyRootOutput();
    void CollimatorEnergyRootOutputWrite(int, int16_t, double);

private:

    TTree *CollimationEnergyTree;

    Int_t id;
    int16_t nucleons;
    Double_t energy;
};

/**
* This class outputs the FLUKA insertion info 
*/
class CollimationFLUKARootOutput
{
public:

    CollimationFLUKARootOutput();

    void FLUKANamesRootOutputWrite(int, char*, int, int);

    std::vector<double> insertion_id;
    std::vector<double> insertion_length;

private:

    TTree *FLUKANamesTree;

    Int_t id;
    Int_t ins_type;
    Char_t name[49];
    Double_t length;

};


#endif

