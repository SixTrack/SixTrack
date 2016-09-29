#ifndef _elasticScatter_h_
#define _elasticScatter_h_

/**
* Include for the vector class for storing the cross section tables.
*/
#include <vector>

/**
* Include for complex numbers
*/
#include <complex>

/**
* Include for the interpolation classes, to interpolate cross section values.
*/
//#include "NumericalUtils/Interpolation.h"
#include "HelperFunctions.h"

namespace ParticleTracking
{

/**
* Class for all things relating to proton-proton elastic scattering.
* This includes, generation of differential cross sections.
* Integration of the differential cross sections.
* Generation of momentum transfer t values for calculation of scattering angles.
* Multiple scattering models
*/
class ppElasticScatter
{
public:
	/**
	* class constructor
	*/
	ppElasticScatter(): Configured(false),Debug(false) {}

//~ppElasticScatter();
	/**
	* Generates the requried differential cross sections and integrates for the specified energy
	*/
	void GenerateTDistribution(double energy);

	/**
	* Generates the elastic differential cross section
	* Places the results into the vectors t and DSig
	* @param energy sqrt s
	*/
	void GenerateDsigDt(double energy);

	/**
	* Integrates the elastic differential cross section
	*/
	void IntegrateDsigDt();

	/**
	* Picks a t value from the generated distribution (including interpolation)
	*/
	double SelectT();

	/**
	* Sets the minimum t value for generation
	* @param tmin the minumum t value to generate
	*/
	void SetTMin(double tmin);

	/**
	* Sets the maximum t value for generation
	*/
	void SetTMax(double);

	/**
	* Sets the step size in t for the differential cross section generation
	*/
	void SetStepSize(double);

	/**
	* Gets the currently set minimum t value
	*/
	double GetTMin() const;

	/**
	* Gets the currently set maximum t value
	*/
	double GetTMax() const;

	/**
	* Gets the currently set step size
	*/
	double GetStepSize() const;

	/**
	* Gets the Integrated elastic cross section
	*/
	double GetElasticCrossSection() const;
	double GetElasticCrossSectionN() const;

	/**
	* Generates the differential cross section at a given t value and energy
	* The energy is the sqrt(s) of the interaction
	*/
	double PomeronScatter(double t_input, double energy, bool em);

	/**
	* Debug toggle - set to true or false to enable/disable debugging output
	*/
	void EnableDebug(bool);

private:
	/**
	* Interpolation classes for the cross section data
	*/
	Interpolation *LinearInterpolation;
	Interpolation *InversionInterpolation;

	/**
	* bool to check if the cross sections have been generated
	*/
	bool Configured;

	/**
	* double to store the minimum t value generate
	*/
	double t_min;

	/**
	* double to store the maximum t value to generate
	*/
	double t_max;

	/**
	* double to store the t step size
	*/
	double step;


	std::vector<double> Uniformt;
	std::vector<double> DSig;
	std::vector<double> IntSig;
	std::vector<double> Sig;

	std::vector<double> DSigN;
//std::vector<double> IntSigN;
	/**
	* The Integrated elastic cross section
	*/
	double SigElastic;
	double SigElasticN;

	/**
	* b slope
	*/
//double b;

	/**
	* Enable Scattering/MC debugging
	*/
	bool Debug;

//double MinSigValue;

	/**
	*
	*	Scattering Functions
	*
	*/


	/**
	*	Proton form factor squared
	*/
	double F1(double t);

	/**
	*	Pomeron trajectory
	*/
	double a0(double t, double *par);

	/**
	*	Pomeron trajectory
	*/
	double a1(double t, double *par);

	/**
	*	 trajectory
	*/
	double a2(double t, double *par);

	/**
	*	 trajectory
	*/
	double a3(double t, double *par);

	/**
	*	Real hard pomeron component
	*/
	double hardpomre(double tnu, double t, double *par);

	/**
	*	Imaginary hard pomeron component
	*/
	double hardpomim(double tnu, double t, double *par);

	double pomre(double tnu, double t, double *par);
	double pomim(double tnu, double t, double *par);

	/**
	*	C=+ve regge trajectory real component
	*/
	double plusre(double tnu, double t, double *par);

	/**
	*	C=+ve regge trajectory imaginary component
	*/
	double plusim(double tnu, double t, double *par);

	/**
	*	C=-ve regge trajectory real component
	*/
	double minusre(double tnu, double t, double *par);

	/**
	*	C=-ve regge trajectory imaginary component
	*/
	double minusim(double tnu, double t, double *par);

	/**
	*	proton-proton two pomeron term
	*/
	std::complex<double> twopom(double tnu, double t, double *par);

	/**
	*	proton-antiproton two pomeron term
	*/
	std::complex<double> twopombar(double tnu, double t, double *par);

	double twopomre(double tnu, double t, double *par);
	double twopomim(double tnu, double t, double *par);

	double twopombarre(double tnu, double t, double *par);
	double twopombarim(double tnu, double t, double *par);

	/**
	*	Triple gluon exchange
	*/
	double ggg(double tt, double *par);

};//End class ppElasticScatter

//typedef std::complex<double> std::complex<double>;

//std::complex<double> Cadd(std::complex<double> a, std::complex<double> b);	//Complex addition
//std::complex<double> Csub(std::complex<double> a, std::complex<double> b);	//Complex subtraction
//std::complex<double> Cmul(std::complex<double> a, std::complex<double> b);	//Complex multiplation
//std::complex<double> Complex_t(double re, double im);				//Constructor
//std::complex<double> Conjg(std::complex<double> z);				//Conjugate
//std::complex<double> Cdiv(std::complex<double> a, std::complex<double> b);	//Complex division
//double Cabs(std::complex<double> z);						//Absolute
//std::complex<double> Csqrt(std::complex<double> z);				//Square root
//std::complex<double> RCmul(double x, std::complex<double> a);			//Real*complex
std::complex<double> RCdiv(double x, std::complex<double> a);			//Real/complex?
//std::complex<double> Cexp(std::complex<double> a);
//std::complex<double> RCadd(double x, std::complex<double> a);

}//End namespace ParticleTracking
#endif
