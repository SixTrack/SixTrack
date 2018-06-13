/**
* Include for the Elastic scattering class
*/
//#include "Collimators/ElasticScatter.h"
#include "ElasticScatter.h"

/**
* Include for the math headers - required for exp, sin, other mathmatical functions
*/
#include <cmath>

/**
* Include for io - std::cout etc
*/
#include <iostream>

/**
* Include for file output
*/
#include <fstream>

/**
* Include for the max() algorithm
*/
#include <algorithm>

/**
* Include for assorted numerial constants
*/
//#include "NumericalUtils/NumericalConstants.h"

/**
* Include for assorted Physical constants
*/
//#include "NumericalUtils/PhysicalConstants.h"

/**
* Include for assorted Physical units
*/
//#include "NumericalUtils/PhysicalUnits.h"

/**
* Include for the random number generator
*/
//#include "Random/RandomNG.h"

/**
* Pulls in <complex> and Complex std::complex<double> typedef
*/
//#include "NumericalUtils/Complex.h"

#include "HelperFunctions.h"

namespace ParticleTracking
{

/**
* Form factor^2 fit to experimental data
*/
double Q[3]= {0.26,0.56,0.18};
double q[3]= {8.38,3.78,1.36};

/**
* trajectory slope parameters.
*/
double apr0=0.1, apr1=0.25, apr2=0.821595, apr3=0.904556;

double par[12];

/**
* Sets the minimum t value for generation
*/
void ppElasticScatter::SetTMin(double tmin)
{
	t_min = tmin;
}

/**
* Sets the maximum t value for generation
*/
void ppElasticScatter::SetTMax(double tmax)
{
	t_max = tmax;
}

/**
* Sets the step size in t for the differential cross section generation
*/
void ppElasticScatter::SetStepSize(double StepSize)
{
	step = StepSize;
}

/**
* Gets the currently set minimum t value
*/
double ppElasticScatter::GetTMin() const
{
	return t_min;
}

/**
* Gets the currently set maximum t value
*/
double ppElasticScatter::GetTMax() const
{
	return t_max;
}

/**
* Gets the currently set step size
*/
double ppElasticScatter::GetStepSize() const
{
	return step;
}

/**
* Gets the Integrated elastic cross section
*/
double ppElasticScatter::GetElasticCrossSection() const
{
	return SigElastic;
}

/**
* Gets the Integrated elastic cross section
*/
double ppElasticScatter::GetElasticCrossSectionN() const
{
	return SigElasticN;
}

/**
* Debug toggle - set to true or false to enable/disable debugging output
*/
void ppElasticScatter::EnableDebug(bool flag)
{
	Debug = flag;
}

/**
* Generates the requried differential cross sections and integrates for the specified energy
*/
void ppElasticScatter::GenerateTDistribution(double energy)
{
	if(!Configured)
	{
		Uniformt = new std::vector<double>;
		DSig = new std::vector<double>;
		DSigN = new std::vector<double>;
		/*
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*******   Generating pp elastic differential cross section   ******************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		*/
		GenerateDsigDt(energy);
		/*
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*************   Integrating differential cross section   **********************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		*/
		IntegrateDsigDt();
		/*
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*************   Configuration generation done!   ******************************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		*/
		Configured = true;
		delete Uniformt;
		delete DSig;
		delete DSigN;
	}
}

ppElasticScatter::~ppElasticScatter()
{
	if(LinearInterpolation)
	{
		delete LinearInterpolation;
	}
}

/**
* Generates the elastic differential cross section
* Places the results into the vectors t and DSig
* @param energy beam energy
*/
void ppElasticScatter::GenerateDsigDt(double energy)
{
	std::cout << "Call Generate DsigDt " << std::endl;

	// Values from James Molson's Thesis

	par[0] = 0.106231;		//eps0
	par[1] = 0.0972043;		//eps1
	par[2] = -0.510662;		//eps2
	par[3] = -0.302082;		//eps3
	par[4] = 228.359;		//X0
	par[5] = 193.811;		//X1
	par[6] = 518.686;			//X2
	par[7] = 10.7843;			//X3
	par[8] = 0.521223;		//lambda
	par[9] = 5.02965;		//t0 for ggg
	par[10] = 0.0449029;		//alpha pomeron
	par[11] = 0.278037;		//alpha hard pomeron


	unsigned int nSteps = (t_max - t_min) / step;
	Uniformt->clear();
	DSig->clear();
	DSigN->clear();
	Uniformt->reserve(nSteps);
	DSig->reserve(nSteps);
	DSigN->reserve(nSteps);

	/**
	* Get the energy of this interaction
	*/
	double s = (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * energy) + (2 * pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV,2));

	//double ppel = 0.007 * pow((energy / 450.0),0.04792);
	//std::cout << "Sixtrack elastic cross section: " << ppel * 1000.0 << " mb" << std::endl;

	//b = 8.5 + 1.086 * log(sqrt(s));
	//std::cout << "sixtrack b slope gradient: " << b << std::endl;

	double sqrts = sqrt(s);
	std::cout << "Using " << nSteps << " bins and sqrt s: " << sqrts << std::endl;

	if(!Debug)
	{
		for(unsigned int n = 0; n <= nSteps; n++)
		{
			Uniformt->push_back((static_cast<double>(n) * step) + t_min);
			DSig->push_back(PomeronScatter((*Uniformt)[n],sqrts,true));
			DSigN->push_back(PomeronScatter((*Uniformt)[n],sqrts,false));
			//Old function
			//DSig.push_back((ppel*b)*exp(-b*Uniformt[n]));
			//std::cout << DSig[n] << "\t" << Uniformt[n] << std::endl;
		}
	}
	else
	{
		std::ofstream *dSigDebug = new std::ofstream("DSigDebugLog");
		for(unsigned int n = 0; n <= nSteps; n++)
		{
			Uniformt->push_back((static_cast<double>(n) * step) + t_min);
			DSig->push_back(PomeronScatter((*Uniformt)[n],sqrts,true));
			DSigN->push_back(PomeronScatter((*Uniformt)[n],sqrts,false));
			(*dSigDebug) << (*Uniformt)[n] << "\t" << (*DSig)[n] << std::endl;

		}
	}
}

/**
* Generates the elastic differential cross section
* Places the results into the vectors t and DSig
*/
void ppElasticScatter::IntegrateDsigDt()
{
	unsigned int nSteps = Uniformt->size();
	std::vector<double> Sig;
	Sig.reserve(nSteps);

	//Add the 0.0 value first!
	std::vector<double> IntSig;
	IntSig.reserve(nSteps);
	IntSig.push_back(0.0);

	SigElastic = 0;
	SigElasticN = 0;
	/**
	* Integrate the elastic differential cross section
	*/
	std::vector<double>::iterator itr;
	itr = DSig->begin()+1;
	std::vector<double>::iterator itrN;
	itrN = DSigN->begin()+1;
	while(itr != DSig->end())
	{
		SigElastic += ((*itr) * step);
		SigElasticN += ((*itrN) * step);
		IntSig.push_back(SigElastic);
		itr++;
		itrN++;
	}

	std::cout << "Elastic Cross section (with peak): " << SigElastic * 1000 << " mb" << std::endl;
	std::cout << "Elastic Cross section (without peak): " << SigElasticN * 1000 << " mb" << std::endl;
	std::cout << "Sixtrack Elastic Cross section: " << 7 * pow((7000 / 450),0.04792) << " mb" << std::endl;

	std::ofstream* ofile;
	std::ofstream* SigmaDistributionFile;
	itr = IntSig.begin()+1;
	if(Debug)
	{
		ofile = new std::ofstream("IntegratedSigNormalized");
		ofile->precision(16);
		(*ofile) << (*Uniformt)[0] << "\t" << IntSig[0] << std::endl;

		//Switch to normalized values to make our life easier
		for(unsigned int n = 1; n < nSteps; n++)
		{
			IntSig[n] /= SigElastic;
			(*ofile) << (*Uniformt)[n] << "\t" << IntSig[n] << std::endl;
		}
	}
	else
	{
		//Switch to normalized values to make our life easier
		while(itr != IntSig.end())
		{
			(*itr) /= SigElastic;
			itr++;
		}
	}

	InversionInterpolation = new Interpolation(IntSig, *Uniformt);

	if(Debug)
	{
		SigmaDistributionFile = new std::ofstream("SigmaTDistribution");
		SigmaDistributionFile->precision(16);
	}

	/*
	* In our normalized distribtion, 0 sigma = t_min
	* 1 sigma = t_max
	* Thus the first entry should be t_min
	*/
	Sig.push_back(t_min);

	if(Debug)
	{
		(*SigmaDistributionFile) << 0.0 << "\t" << t_min << std::endl;
	}

	double sig_gen = 0;

	for(unsigned int n=1; n <nSteps; n++)
	{
		double target = (static_cast<double>(n) / nSteps);
		try
		{
			sig_gen = (*InversionInterpolation)(target);
		}
		catch(Interpolation::BadRange& error)
		{
			std::cerr << "Bad Range in interpolation - requested: " << error.value << " but valid range is from " << error.valid_range.lower << " to "  << error.valid_range.upper << std::endl;
			std::cerr << "error in entry: " << n << " with total " << nSteps << std::endl;
			throw;
		}

		if(Debug)
		{
			(*SigmaDistributionFile) << target << "\t" << sig_gen << std::endl;
		}

		Sig.push_back(sig_gen);
	}

	/*
	* In our normalized distribtion, 0 sigma = t_min
	* 1 sigma = t_max
	* Thus the last entry should be t_max
	*/
	Sig.push_back(t_max);

	if(Debug)
	{
		(*SigmaDistributionFile) << 1.0 << "\t" << t_max << std::endl;
	}

	LinearInterpolation = new Interpolation(Sig, 0, (1.0/nSteps));    // Interpolation of equally spaced data points

	if(Debug)
	{
		SigmaDistributionFile->close();
		ofile->close();

		//Free up memory
		delete SigmaDistributionFile;
		delete ofile;
	}

	delete InversionInterpolation;
	//delete LinearInterpolation;
}

/**
* Picks a t value from the generated distribution (including interpolation)
*/
double ppElasticScatter::SelectT()
{
	double SigValue = RandomNG::uniform (0, 1.0);
	double t = (*LinearInterpolation)(SigValue);
	/*
		if(Debug)
		{
			std::ofstream *out = new std::ofstream("GeneratedTValues",std::ios_base::app);
			(*out) << SigValue << "\t" << t << std::endl;
		}
	*/
	return t;
}

double ppElasticScatter::PomeronScatter(double tt, double sqrt_s, bool em)
{
	//std::cout << "Call PomeronScatter elastic " << std::endl;
	apr1=par[10];
	apr0=par[11];
	double mu = 0.93827203;
	double m = 0.93827203;

	double real, imag;//t=0
	double tnu=sqrt_s*sqrt_s-m*m-mu*mu-0.5*tt;

	double bs1 = 8.1;
	double bs2 = 1.2;
	bs1 = 8.5;
	bs2 = 1.086;

	//Em coupling ~ 1/137
	double alpha = 0.0072973525376;
	if(!em)
	{
		alpha=0;
	}

	//Euler constant
	double gamma = 0.577215664901532861;
	double bslope = bs1 + bs2*log(sqrt_s);
	double ppC;
	//double phi = -(gamma + log(0.5*bslope));
	double phi = -(gamma + log(tt*bslope/2) + log(1 +  8/(bslope*0.71)) + ((4*tt/0.71) * log((4*tt)/0.71)) + (2*tt/0.71));
	real=(hardpomre(tnu,tt,par) + pomre(tnu,tt,par) + plusre(tnu,tt,par) + minusre(tnu,tt,par) + par[8]*twopomre(tnu,tt,par)) / tnu + ggg(tt,par);
	imag=(hardpomim(tnu,tt,par) + pomim(tnu,tt,par) + plusim(tnu,tt,par) + minusim(tnu,tt,par) + par[8]*twopomim(tnu,tt,par)) / tnu;

	double ac = -alpha * F1(tt)/tt;
	real *=sqrt(0.25*0.389);
	imag *=sqrt(0.25*0.389);
	double rho = real/imag;
	//ac = 0;
	ppC = (4.0 * pow(pi,2.0) * pow(ac,2.0) + pow(real,2.0) + pow(imag,2.0) + 2.0*(rho+alpha*phi) * 2.0*pi*ac*imag) / (4.0*pi);

	return (ppC/1000.0);	//Convert to barns
}

/**
*
*	Scattering Function internals
*
*/

/**
*	Proton form factor squared
*/
double ppElasticScatter::F1(double tt)
{
	double S=0;
	for(int i=0; i<3; i++)
	{
		S += Q[i] * exp(-q[i]*tt);
	}
	return S;
}

/**
*	 trajectory
*/
double ppElasticScatter::a0(double tt, double *par)
{
	return (1.0+par[0]-apr0*tt);
}

/**
*	 trajectory
*/
double ppElasticScatter::a1(double tt, double *par)
{
	return 1.0+par[1]-apr1*tt;
}

/**
*	 trajectory
*/
double ppElasticScatter::a2(double tt, double *par)
{
	return 1.0+par[2]-apr2*tt;
}

/**
*	 trajectory
*/
double ppElasticScatter::a3(double tt, double *par)
{
	return 1.0+par[3]-apr3*tt;
}
//______________________________________________________________________________
double ppElasticScatter::hardpomre(double tnu, double tt, double *par)
{
	return -par[4]*F1(tt)*cos(0.5*pi*a0(tt,par))*pow(tnu*apr0,a0(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::hardpomim(double tnu, double tt, double *par)
{
	return par[4]*F1(tt)*sin(0.5*pi*a0(tt,par))*pow(tnu*apr0,a0(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::pomre(double tnu, double tt, double *par)
{
	return -par[5]*F1(tt)*cos(0.5*pi*a1(tt,par))*pow(tnu*apr1,a1(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::pomim(double tnu, double tt, double *par)
{
	return par[5]*F1(tt)*sin(0.5*pi*a1(tt,par))*pow(tnu*apr1,a1(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::plusre(double tnu, double tt, double *par)
{
	return -par[6]*F1(tt)*cos(0.5*pi*a2(tt,par))*pow(tnu*apr2,a2(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::plusim(double tnu, double tt, double *par)
{
	return par[6]*F1(tt)*sin(0.5*pi*a2(tt,par))*pow(tnu*apr2,a2(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::minusre(double tnu, double tt, double *par)
{
	return -par[7]*F1(tt)*sin(0.5*pi*a3(tt,par))*pow(tnu*apr3,a3(tt,par));
}
//______________________________________________________________________________
double ppElasticScatter::minusim(double tnu, double tt, double *par)
{
	return -par[7]*F1(tt)*cos(0.5*pi*a3(tt,par))*pow(tnu*apr3,a3(tt,par));
}
//______________________________________________________________________________
std::complex<double> ppElasticScatter::twopom(double tnu, double tt, double *par)
{
	int i,j, u,v;
	double T;
	double X[4]= {par[4],par[5],par[6],par[7]};
	double apr[4]= {apr0,apr1,apr2,apr3};
	double e[4]= {par[0],par[1],par[2],par[3]};

	std::complex<double> D[4][3],L[4],LL[4];
	std::complex<double> S(0,0),U,V,W,Z;

	for(u=0; u<4; u++)
	{
		L[u]=std::complex<double>(apr[u]*log(apr[u]*tnu),-0.5*apr[u]*pi);
		LL[u]=((1.0+e[u])/apr[u])*L[u];
		for(i=0; i<3; i++)
		{
			D[u][i]=q[i]+L[u];
		}
	}
	for(u=0; u<4; u++)
	{
		for(v=0; v<4; v++)
		{
			for(i=0; i<3; i++)
			{
				for(j=0; j<3; j++)
				{
					U=D[u][i]*D[v][j];
					V=D[u][i]+D[v][j];
					//W=Cdiv(U,V);
					W=U/V;
					W*=-tt;
					W=exp(W);
					Z=RCdiv(Q[i]*Q[j],V);
					Z*=W;
					T=X[u]*X[v]/(16.0*pi*tnu);
					Z*=T;
					V=LL[u]+LL[v];
					V=exp(V);
					Z*=V;
					if(u==3)
					{
						Z*=Complex(0,1.0);
					}
					if(v==3)
					{
						Z*=Complex(0,1.0);
					}
					S+=Z;
				}
			}
		}
	}
	S*=Complex(0,1.0);
	return S;
}
//______________________________________________________________________________
std::complex<double> ppElasticScatter::twopombar(double tnu, double tt, double *par)
{
	int i,j, u,v;
	double T;
	double X[4]= {par[4],par[5],par[6],-par[7]};
	double apr[4]= {apr0,apr1,apr2,apr3};
	double e[4]= {par[0],par[1],par[2],par[3]};

	std::complex<double> D[4][3],L[4],LL[4];
	std::complex<double> S(0,0),U,V,W,Z;

	for(u=0; u<4; u++)
	{
		L[u]=std::complex<double>(apr[u]*log(apr[u]*tnu),-0.5*apr[u]*pi);
		LL[u]=((1.0+e[u])/apr[u])*L[u];
		for(i=0; i<3; i++)
		{
			D[u][i]=q[i]+L[u];
		}
	}
	for(u=0; u<4; u++)
	{
		for(v=0; v<4; v++)
		{
			for(i=0; i<3; i++)
			{
				for(j=0; j<3; j++)
				{
					U=D[u][i]*D[v][j];
					V=D[u][i]+D[v][j];
					//W=Cdiv(U,V);
					W=U/V;
					W*=-tt;
					W=exp(W);
					Z=RCdiv(Q[i]*Q[j],V);
					Z*=W;
					T=X[u]*X[v]/(16.0*pi*tnu);
					Z*=T;
					V=LL[u]+LL[v];
					V=exp(V);
					Z*=V;
					if(u==3)
					{
						Z*=Complex(0,1.0);
					}
					if(v==3)
					{
						Z*=Complex(0,1.0);
					}
					S*=Z;
				}
			}
		}
	}
	S*=Complex(0,1.0);
	return S;
}
//______________________________________________________________________________
double ppElasticScatter::twopomre(double tnu, double tt, double *par)
{
	//std::complex<double> res = twopom(tnu,tt,par);
	//return res.real();
	return twopom(tnu,tt,par).real();
}
//______________________________________________________________________________
double ppElasticScatter::twopomim(double tnu, double tt, double *par)
{
//	std::complex<double> res = twopom(tnu,tt,par);
//	return res.imag();
	return twopom(tnu,tt,par).imag();
}
//______________________________________________________________________________
double ppElasticScatter::twopombarre(double tnu, double tt, double *par)
{
//	std::complex<double> res = twopombar(tnu,tt,par);
//	return res.real();
	return twopombar(tnu,tt,par).real();
}
//______________________________________________________________________________
double ppElasticScatter::twopombarim(double tnu, double tt, double *par)
{
	//std::complex<double> res = twopombar(tnu,tt,par);
	//return res.imag();
	return twopombar(tnu,tt,par).imag();
}

/**
*	Triple gluon exchange
*/
double ppElasticScatter::ggg(double tt, double *par)
{
	double CCC=sqrt(0.09*16.0*pi/0.389);
	if(tt>par[9])
	{
		return CCC/pow(tt,4.0);
	}
	else
	{
		return CCC*exp(4.0-4.0*tt/par[9])/pow(par[9],4.0);
	}
}

// complex.c
std::complex<double> RCdiv(double x, std::complex<double> a)
{
	double den=(a.real())*(a.real())+(a.imag())*(a.imag());
	std::complex<double> c(x*(a.real())/den,-x*(a.imag())/den);
	return c;
}

/* (C) Copr. 1986-92 Numerical Recipes Software 06,. */

}//End namespace ParticleTracking
