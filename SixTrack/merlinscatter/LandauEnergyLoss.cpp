#include "HelperFunctions.h"
#include "LandauEnergyLoss.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

//(mat,p,rlen,m_dpodx)
//Advanced energy loss
//void ScatteringModel::EnergyLoss(PSvector& p, double x, Material* mat, double E0)
/**
* We need to know the nucleus we are working with here in "Merlin".
The parameters must be sent over from sixtrack.
For the electron density (and thus the plasma energy), we need the at
*/
double EnergyLoss(double p, double ElectronDensity, double PlasmaEnergy, double MeanIonisationEnergy, double StepSize)
{
	//Convert I to eV (comes in as GeV)
	double I = MeanIonisationEnergy/eV;

	//Check p units (should be GeV)
	double gamma = p/(ProtonMassMeV*MeV);
	double beta = sqrt(1 - ( 1 / (gamma*gamma)));


	double RandomLandau = RandomNG::landau();

	double tmax = (2*ElectronMassMeV * beta * beta * gamma * gamma ) / (1 + (2 * gamma * (ElectronMassMeV/ProtonMassMeV)) + pow((ElectronMassMeV/ProtonMassMeV),2))*MeV;

	static const double xi1 = 2.0 * pi * pow(ElectronRadius,2) * ElectronMass * pow(SpeedOfLight,2);
	double xi0 = xi1 * ElectronDensity;
	double xi = (xi0 * StepSize /(beta*beta)) / ElectronCharge * (eV/MeV);

	double C = 1 + 2*log(I/(PlasmaEnergy/eV));
	double C1 = 0;
	double C0 = 0;

	if((I/eV) < 100)
	{
		if(C <= 3.681)
		{
			C0 = 0.2;
			C1 = 2.0;
		}
		else
		{
			C0 = 0.326*C - 1.0;
			C1 = 2.0;
		}
	}
	else	//I >= 100eV
	{
		if(C <= 5.215)
		{
			C0 = 0.2;
			C1 = 3.0;
		}
		else
		{
			C0 = 0.326*C - 1.5;
			C1 = 3.0;
		}
	}
	double delta = 0;

	//Density correction
	double ddx = log10(beta*gamma);
	if(ddx > C1)
	{
		delta = 4.606*ddx - C;
	}
	else if(ddx >= C0 && ddx <= C1)
	{
		double m = 3.0;
		double xa = C /4.606;
		double a = 4.606 * (xa - C0) / pow((C1-C0),m);
		delta = 4.606*ddx -C + a*pow((C1 - ddx),m);
	}
	else
	{
		delta = 0.0;
	}

//	double tcut = 2.0*MeV;
//	tcut = tmax;

	//Mott Correction
	double G = pi*FineStructureConstant*beta/2.0;
	double q = (2*(tmax/MeV)*(ElectronMassMeV) )/(pow((0.843/MeV),2));
	double S = log(1+q);
	double L1 = 0.0;
	double yL2 = FineStructureConstant/beta;

	double L2sum = 1.202001688211;	//Sequence limit calculated with mathematica for beta = 1
	double L2 = -yL2*yL2*L2sum;

	double F = G - S + 2*(L1 + L2);
	//double deltaE = xi * (log(2 * ElectronMassMeV * beta*beta * gamma*gamma * (tcut/MeV)/pow(I/MeV,2)) - (beta*beta)*(1 + ((tcut/MeV)/(tmax/MeV))) - delta + F - 1.0 - euler);
	double deltaE = xi * (log(2 * ElectronMassMeV * beta*beta * gamma*gamma * xi /pow(I/MeV,2)) - (beta*beta) - delta + F + 0.20);

	double dp = ((xi * RandomLandau) - deltaE) * MeV;

	return dp;
}

/*
double Material::CalculateElectronDensity()
{
	return AtomicNumber * Avogadro * Density * 0.001 / (AtomicMass * pow(centimeter,3)); // n_e m^-3 (1e6 conversion from cm^3)
}
double Material::CalculatePlasmaEnergy()
{
	return (PlanckConstantBar * sqrt((ElectronDensity * pow(ElectronCharge,2)) / (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
}

double Material::CalculateMeanExcitationEnergy()
{
	return AtomicNumber * 10.0 *eV;
}
*/
