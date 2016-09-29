#include <iostream>
#include "DiffractiveScatter.h"
#include "ElasticScatter.h"
#include "RandomNG/RandomNG.h"

using namespace ParticleTracking;

bool MScatterConfigured = false;
double p_store;
ppDiffractiveScatter* ppDiffractiveScatter_ptr = nullptr;
ppElasticScatter* ppElasticScatter_ptr = nullptr;

/**
* Configuration function for the Merlin scattering physics.
* Normally called once at the start of the tracking run.
* This will configure all the sub classes required.
* @param[in] Plab The reference momentum of the proton beam.
*/
extern "C" void merlinscatter_setup_(double* Plab, int* seed)
{
	//std::cout << "merlinscatter_setup(), plab = "<< *Plab << std::endl;

	//Initialise Random number generator
	if(!MScatterConfigured)
	{
		RandomNG::init(*seed);

		//Store the beam momentum so we have access to it later.
		p_store = *Plab;

		//First configure the diffractive

		//If we already have a class created, delete it and reset
		//if(!ppDiffractiveScatter_ptr)
		//{
		std::cout << std::endl << "merlinscatter_setup(), configure Single Diffractive scattering start"<< std::endl;

		double Mproton = 0.938272013;
		double Mpion = 0.1349766;
		double s = (2*pow(Mproton,2)+2*Mproton*p_store);
		double Mmin2 = pow(Mproton+Mpion,2);
		double xi_th = Mmin2/s; // (M_p + M_pion)^2/s

		ppDiffractiveScatter_ptr = new ppDiffractiveScatter();

		ppDiffractiveScatter_ptr->SetTMin(0.0001);
		ppDiffractiveScatter_ptr->SetTMax(4);
		ppDiffractiveScatter_ptr->SetTStepSize(1e-4);

		ppDiffractiveScatter_ptr->SetXiMin(xi_th);//Threshould at (M_proton + M_pion)^2/s
		ppDiffractiveScatter_ptr->SetXiMax(0.12);
		ppDiffractiveScatter_ptr->SetXiStepSize(1e-6);

		ppDiffractiveScatter_ptr->GenerateDistribution(s);

		std::cout << "merlinscatter_setup(), configure Single Diffractive scattering end" << std::endl << std::endl;
		//}

		//Now the elastic
		//if(!ppElasticScatter_ptr)
		//{
		std::cout << "merlinscatter_setup(), configure Elastic scattering start" << std::endl;

		ppElasticScatter_ptr = new ppElasticScatter();
		ppElasticScatter_ptr->SetTMin(1e-4);
		ppElasticScatter_ptr->SetTMax(1.0);
		ppElasticScatter_ptr->SetStepSize(1e-4);
		ppElasticScatter_ptr->GenerateTDistribution(p_store);

		std::cout << "merlinscatter_setup(), configure Elastic scattering end" << std::endl << std::endl;

		MScatterConfigured = true;
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
	//std::cout << "merlinscatter_setdata()" << std::endl;

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
	//std::cout << "End merlinscatter_setdata()" << std::endl;
}

/**
* Selects a t value for elastic scattering
*/
extern "C" void merlinscatter_get_elastic_t_(double* gettran)
{
	*gettran = ppElasticScatter_ptr->SelectT();
}

/**
* Selects a t value for single diffractive scattering
*/
extern "C" void merlinscatter_get_sd_t_(double* gettran)
{
	*gettran = ppDiffractiveScatter_ptr->Select().first;
}

/**
* Selects a xi value for single diffractive scattering
*/
extern "C" void merlinscatter_get_sd_xi_(double* xm2)
{
	std::pair<double,double>TM = ppDiffractiveScatter_ptr->Select();
	double m_rec = TM.second;
	*xm2 = m_rec;
//	double com_sqd = (2 * ProtonMassMeV * MeV * p_store) + (2 * ProtonMassMeV * MeV * ProtonMassMeV * MeV);
//	double dp = m_rec * m_rec * E / com_sqd;
}
