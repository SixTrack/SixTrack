//For scattering information output
#include <iostream>

//For single diffraction scattering
#include "DiffractiveScatter.h"

//For elastic scattering
#include "ElasticScatter.h"

//For ionization energy losses with landau straggling
#include "LandauEnergyLoss.h"

//For the Merlin random number generator
#include "RandomNG/RandomNG.h"

using namespace ParticleTracking;

//A bool to check if the scattering classes are configured yet
bool MScatterConfigured = false;

//An internal store of the reference momentum
double p_store;

//The pointer to the diffractive scattering class
ppDiffractiveScatter* ppDiffractiveScatter_ptr = nullptr;

//The pointer to the elastic scattering class
ppElasticScatter* ppElasticScatter_ptr = nullptr;

/**
* Configuration function for the Merlin scattering physics.
* Normally called once at the start of the tracking run.
* This will configure all the sub classes required.
* @param[in] Plab The reference momentum of the proton beam.
* @param[in] seed The seed value to initialise the RNG with.
*/
extern "C" void merlinscatter_setup_(double* Plab, int* seed)
{
	if(!MScatterConfigured)
	{
		std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
		std::cout << "Configuring Merlin RNG with the seed " << *seed << std::endl;

		//Initialise Random number generator
		RandomNG::init(*seed);

		//Store the beam momentum so we have access to it later.
		p_store = *Plab;

		std::cout << "Using input beam energy of " << *Plab << " GeV" << std::endl;

		//First configure the diffractive scattering
		std::cout << std::endl << "merlinscatter_setup(), configure Single Diffractive scattering start"<< std::endl;

		double Mproton = PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV;
		double Mpion = 0.1349766;
//		double s = (2*pow(Mproton,2)+2*Mproton*p_store);
		double s = (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * *Plab) + (2 * pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV,2));
		double Mmin2 = pow(Mproton+Mpion,2);
		double xi_th = Mmin2/s; // (M_p + M_pion)^2/s

		//Make our diffractive scattering class
		ppDiffractiveScatter_ptr = new ppDiffractiveScatter();

		//Configure the t values and step size to generate (can be adjusted as needed).
		ppDiffractiveScatter_ptr->SetTMin(0.0001);
		ppDiffractiveScatter_ptr->SetTMax(4);
		ppDiffractiveScatter_ptr->SetTStepSize(1e-4);

		//Configure the xi values and step size to generate (can be adjusted as needed).
		ppDiffractiveScatter_ptr->SetXiMin(xi_th);//Threshould at (M_proton + M_pion)^2/s
		ppDiffractiveScatter_ptr->SetXiMax(0.12);
		ppDiffractiveScatter_ptr->SetXiStepSize(1e-6);

		//Make the differential cross section look up table
		ppDiffractiveScatter_ptr->GenerateDistribution(*Plab);

		std::cout << "merlinscatter_setup(), configure Single Diffractive scattering end" << std::endl << std::endl;

		//Now configure the elastic scattering.
		std::cout << "merlinscatter_setup(), configure Elastic scattering start" << std::endl;

		//Make our elastic scattering class
		ppElasticScatter_ptr = new ppElasticScatter();

		//Configure the t values and step size to generate (can be adjusted as needed).
		ppElasticScatter_ptr->SetTMin(1e-4);
		ppElasticScatter_ptr->SetTMax(1.0);
		ppElasticScatter_ptr->SetStepSize(1e-4);

		//Make the differential cross section look up table
		ppElasticScatter_ptr->GenerateTDistribution(*Plab);

		std::cout << "merlinscatter_setup(), configure Elastic scattering end" << std::endl;

		//Set our bool to true so that we do not perform the expensive cross section generation again
		MScatterConfigured = true;

		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	}
}

/**
* Calculates the total cross sections for each process and returns them to the sixtrack configuration routine
* @param[out] pptot The total p-p cross section
* @param[out] ppel The total p-p elastic cross section
* @param[out] ppsd The total p-p single diffractive cross section
*/
extern "C" void merlinscatter_setdata_(double* pptot, double* ppel, double* ppsd)
{
	double Mproton = 0.938272013;
	double s = (2*pow(Mproton,2)+2*Mproton*p_store);

	//Parameters for the Merlin total pp cross section, just taken from the current PDG equation.
	const double Z_pp = 35.4548e-3;
	const double B_pp = 0.30810e-3;
	const double Y1_pp = 42.5323e-3;
	const double Y2_pp = 33.3433e-3;
	const double eta1 = 0.45817;
	const double eta2 = 0.545;
	const double s0 = 28.998225; //*PhysicalUnits::GeV*PhysicalUnits::GeV;
	const double s1 =1; //*PhysicalUnits::GeV*PhysicalUnits::GeV;

	//The Merlin/PDG total pp cross section.
	double sig_tot_m = Z_pp + B_pp*pow(log (s/s0),2.0) + Y1_pp * pow(s1/s,eta1) -Y2_pp * pow(s1/s,eta2);
	//std::cout << "Setting total cross section: " << sig_tot_m << std::endl;

	//The total pp elastic cross section at the configured energy.
	double sig_el_m = ppElasticScatter_ptr->GetElasticCrossSection();
	//std::cout << "Setting elastic total cross section: " << sig_el_m << std::endl;

	//The total pp single diffraction cross section at the configured energy.
	double sig_sd_m = ppDiffractiveScatter_ptr->GetDiffractiveCrossSection();
	//std::cout << "Setting single diffractive total cross section: " << sig_sd_m << std::endl;

	*pptot = sig_tot_m;
	*ppel  = sig_el_m;
	*ppsd  = sig_sd_m;
}

/**
* Selects a t value for elastic scattering
* @param[out] gettran The generated momentum transfer value.
*/
extern "C" void merlinscatter_get_elastic_t_(double* gettran)
{
	*gettran = ppElasticScatter_ptr->SelectT();
}

/**
* Selects a t value for single diffractive scattering
* @param[out] gettran The generated momentum transfer value.
*/
extern "C" void merlinscatter_get_sd_t_(double* gettran)
{
	*gettran = ppDiffractiveScatter_ptr->Select().first;
}

/**
* Selects a xi value for single diffractive scattering
* @param[out] xm2 The generated xi value.
*/
extern "C" void merlinscatter_get_sd_xi_(double* xm2)
{
	std::pair<double,double>TM = ppDiffractiveScatter_ptr->Select();
	double m_rec = TM.second;
	*xm2 = m_rec;
//	double com_sqd = (2 * ProtonMassMeV * MeV * p_store) + (2 * ProtonMassMeV * MeV * ProtonMassMeV * MeV);
//	double dp = m_rec * m_rec * E / com_sqd;
}

extern "C" void merlinscatter_calc_ion_loss_(double* p, double* ElectronDensity, double* PlasmaEnergy, double* MeanIonisationEnergy, double* StepSize, double* LostE)
{
/*
	std::cout << "Energy: " << *p << std::endl;
	std::cout << "ElectronDensity: " << *ElectronDensity << std::endl;
	std::cout << "PlasmaEnergy: " << *PlasmaEnergy << std::endl;
	std::cout << "MeanIonisationEnergy: " << *MeanIonisationEnergy << std::endl;
	std::cout << "StepSize: " << *StepSize << std::endl;
*/
	*LostE = EnergyLoss(*p, *ElectronDensity, *PlasmaEnergy, *MeanIonisationEnergy, *StepSize);
}

