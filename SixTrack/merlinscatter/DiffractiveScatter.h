#ifndef _DiffractiveScatter_h_
#define _DiffractiveScatter_h_

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
* Class for all things relating to proton-proton single diffractive scattering.
* This includes, generation of differential cross sections.
* Integration of the differential cross sections.
* Generation of momentum transfer t values for calculation of scattering angles.
* Generation of mass loss values for calculation of momentum changes.
* Multiple scattering models
*/
class ppDiffractiveScatter
{
public:
	/**
	* class constructor
	*/
	ppDiffractiveScatter(): Configured(false),Debug(false) {}
	~ppDiffractiveScatter();
	/**
	* Generates the requried differential cross sections and integrates for the specified energy
	*/
	void GenerateDistribution(double energy);

	/**
	* Generates the elastic differential cross section
	* Places the results into the vectors t and DSig
	* @param energy sqrt s
	*/
	void GenerateDsigDtDxi(double energy);

	/**
	* Integrates the elastic differential cross section
	*/
	void IntegrateDsigDtDxi();

	/**
	* Picks a t value from the generated distribution (including interpolation)
	*/
	double SelectT();

	/**
	* Picks an xi value from the generated distribution (including interpolation)
	*/
	double SelectXi();

	/**
	* Sets the minimum t value for generation
	* @param tmin the minumum t value to generate
	*/
	void SetTMin(double tmin);

	/**
	* Sets the maximum t value for generation
	*/
	void SetTMax(double tmax);

	/**
	* Sets the maximum t value for generation
	* @param tmin the minumum xi value to generate
	*/
	void SetXiMax(double ximax);

	/**
	* Sets the minimum xi value for generation
	* @param tmin the minumum xi value to generate
	*/
	void SetXiMin(double ximin);

	/**
	* Sets the step size in t for the differential cross section generation
	* @param step The step size to generate
	*/
	void SetTStepSize(double step);

	/**
	* Sets the step size in Xi for the differential cross section generation
	* @param step The step size to generate
	*/
	void SetXiStepSize(double step);

	/**
	* Gets the currently set minimum t value
	*/
	double GetTMin() const;

	/**
	* Gets the currently set maximum t value
	*/
	double GetTMax() const;

	/**
	* Gets the currently set minimum xi value
	*/
	double GetXiMin() const;

	/**
	* Gets the currently set maximum xo value
	*/
	double GetXiMax() const;

	/**
	* Gets the currently set t step size
	*/
	double GetTStepSize() const;

	/**
	* Gets the currently set xi step size
	*/
	double GetXiStepSize() const;

	/**
	* Gets the Integrated Single diffractive cross section
	*/
	double GetDiffractiveCrossSection() const;


	/**
	* Debug toggle - set to true or false to enable/disable debugging output
	* @param debug Turn on and off debugging output
	*/
	void EnableDebug(bool debug);

	//Get our scattering (roger)
	std::pair<double,double> Select();

private:
	/**
	* Generates the differential cross section at a given t value and energy
	* The energy is the sqrt(s) of the interaction
	*/
	inline double PomeronScatter(const double t_input,const double xi_input,const double energy) const;
	double PomeronScatter2(double t_input,double xi_input,double energy);

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
	* double to store the minimum t value to generate
	*/
	double t_min;

	/**
	* double to store the maximum t value to generate
	*/
	double t_max;

	/**
	* double to store the minimum xi value to generate
	*/
	double xi_min;

	/**
	* double to store the maximum xi value to generate
	*/
	double xi_max;

	/**
	* double to store the t step size
	*/
	double t_step;

	/**
	* double to store the xi step size
	*/
	double xi_step;

	std::vector<double> UniformT;
	std::vector<double> UniformXi;
	std::vector<double> DSig;
	std::vector<double> IntSig;
	std::vector<double> Sig;

	/**
	* The Integrated Single diffractive cross section
	*/
	double SigDiffractive;

	/**
	* Enable Scattering/MC debugging
	*/
	bool Debug;

	//Roger stuff - array size
	static const int N = 500;

	double xarray[N];
	double tarray[N];

	//s of the interaction
	double ss;

	/**
	*
	*	Scattering Functions
	*
	*/

};//End class ppDiffractiveScatter

}//End namespace ParticleTracking
#endif
