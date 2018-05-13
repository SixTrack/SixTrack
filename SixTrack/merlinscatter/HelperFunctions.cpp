//#include "NumericalUtils/Interpolation.h"
#include "HelperFunctions.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <cassert>


using namespace std;

// Method used for arbitrarily space data
//
class ArbSpacedData : public Interpolation::Method
{
public:

	struct Data
	{
		double x,y;
		Data(double x1=0, double y1=0) : x(x1),y(y1) {}
	};

	ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals);
	virtual double ValueAt(double) const;

private:
	vector<Data> itsData;
};

// Method used for equally space data
//
class EqualSpacedData : public Interpolation::Method
{
public:

	EqualSpacedData(const vector<double> yv, double xm, double delta)
		: yvals(yv),xmin(xm),xmax(xm+(yv.size()-1)*delta),dx(delta)
	{
		assert(dx>0);
	}

	double ValueAt(double x) const;

private:
	vector<double> yvals;
	double xmin;
	double xmax;
	double dx;
};

// Function  used by sort() for ArbSpacedData::Data
//
inline bool sort_x(const ArbSpacedData::Data& d1, const ArbSpacedData::Data& d2)
{
	return d1.x<d2.x;
}

// Class ArbSpacedData implementation
//
ArbSpacedData::ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals)
	: itsData()
{
	assert(xvals.size()==yvals.size());
	itsData.reserve(xvals.size());
	for(size_t i=0; i<xvals.size(); i++)
	{
		itsData.push_back(Data(xvals[i],yvals[i]));
	}
	sort(itsData.begin(),itsData.end(),sort_x);
}

double ArbSpacedData::ValueAt(double x) const
{
	if(x<itsData.front().x || x>itsData.back().x)
	{
		throw Interpolation::BadRange(x,FloatRange(itsData.front().x,itsData.back().x));
	}

	// locate segment by binary search
	size_t ju=itsData.size(),jl=0,jm;

	while((ju-jl)>1)
	{
		jm = (ju+jl)>>1;
		if(x>itsData[jm].x)
		{
			jl=jm;
		}
		else
		{
			ju=jm;
		}
	}

	// use linear interpolation to return value
	double m = (itsData[jl+1].y-itsData[jl].y)/(itsData[jl+1].x-itsData[jl].x);
	return itsData[jl].y+(x-itsData[jl].x)*m;
}

// Class EqualSpacedData implementation
//

double EqualSpacedData::ValueAt(double x) const
{
	// note that we use extrapolation here if x is out of range
	size_t n;
	if(x<xmin)
	{
		n=0;
	}
	else if(x>xmax)
	{
		n=yvals.size()-2;
	}
	else
	{
		n = (x-xmin)/dx;
	}

	double m = (yvals[n+1]-yvals[n])/dx;
	double x0 = xmin+dx*n;
	return yvals[n]+m*(x-x0);
}



// exception
Interpolation::BadRange::BadRange(double x, const FloatRange& r)
	: MerlinException(""),value(x),valid_range(r)
{
	ostringstream buf;
	buf<<value<<" not in interpolation range ("<<r.lower<<","<<r.upper<<')';
	SetMsg(buf.str());
}

Interpolation::Interpolation(const vector<double>& yvals, double xmin, double dx)
	: itsMethod(new EqualSpacedData(yvals,xmin,dx))
{}

Interpolation::Interpolation(const std::vector<double>& xv, const vector<double>& yv)
	: itsMethod(new ArbSpacedData(xv,yv))
{}

Interpolation::~Interpolation()
{
	if(itsMethod)
	{
		delete itsMethod;
	}
}

// Class Utility PhysicalConstants
namespace PhysicalConstants
{
using namespace PhysicalUnits;

//Updated with PDG 2008 values

const double Avogadro = 6.02214179e23;

const double AtomicMassUnit = 931.494028*MeV; //MeV given by PDG 2008

const double SpeedOfLight = 2.99792458e+08*meter/second;

const double ElectronMass = 9.10938215e-31; // kilogram (old = 9.10956e-31)

const double ProtonMass = 1.672621637e-27; // kilogram  (old = 1.67261e-27)

const double PlanckConstant = 6.62606896e-34; // Joule second (old = 6.62620e-34)

const double PlanckConstantBar = PlanckConstant/(2*pi);

const double ElectronCharge = 1.602176487e-19; // Coulomb (old = 1.60219e-19)

const double ElectronMassMeV = 0.510998910; // MeV (old = 0.511005)

const double ProtonMassMeV = 938.272046; // MeV (old = 938.2578)

const double ProtonMassGeV = 0.938272046;

const double FreeSpacePermeability = 16.0e-7*atan(1.0); // Henry per meter

const double FreeSpacePermittivity = 1.0/FreeSpacePermeability/SpeedOfLight/SpeedOfLight; // Farad per meter

const double ElectronRadius = ElectronCharge*ElectronCharge/16.0/atan(1.0)/FreeSpacePermittivity/ElectronMass/SpeedOfLight/SpeedOfLight; // meter

const double ElectronGe = 0.001159652181; //Magnetic moment anomaly e cm (old = 0.00115965219)

const double MuonMassMeV = 105.65836;

const double MuonMass = ProtonMass * MuonMassMeV/ProtonMassMeV; //AtomicMassUnit * 0.11342892;

const double MuonLifetime = 2.1970e-6 * second;	//In the rest frame

const double FineStructureConstant = (ElectronCharge * ElectronCharge * SpeedOfLight * FreeSpacePermeability) / (2 * PlanckConstant);

const double PionZeroMassMeV = 134.976;

double LorentzBeta (double gamma)
{
	return sqrt( 1-(1/(gamma*gamma)) );
}
double LorentzGamma (double beta)
{
	return 1 / (sqrt(1-(beta*beta)));
}
double LorentzGamma (double momentum, double mass)
{
	return sqrt(1 + pow(momentum/(mass*SpeedOfLight),2) );
}
}

namespace PhysicalUnits
{


//Area
const double barn = 1.0e-28;	// m^2

// length
const double meter=1.0;
const double centimeter =1.0e-02*meter;
const double millimeter =1.0e-03*meter;
const double micrometer =1.0e-06*meter;
const double nanometer  =1.0e-09*meter;

// time
const double second = 1;
const double millisecond = 1.0e-03*second;
const double microsecond = 1.0e-06*second;
const double nanosecond = 1.0e-09*second;
const double picosecond = 1.0e-12*second;
const double minute = 60*second;
const double hour = 60*minute;
const double day = 24*hour;
const double year = 365*day;

// energy
const double eV = 1.0e-9;
const double keV = 1.0e+3*eV;
const double MeV = 1.0e+6*eV;
const double GeV = 1.0e+9*eV;
const double TeV = 1.0e+12*eV;

// voltage
const double Volt = 1.0e-9;
const double kV = 1.0e+3*Volt;
const double MV = 1.0e+6*Volt;
const double GV = 1.0e+9*Volt;
const double TV = 1.0e+12*Volt;

// Frequency
const double Hz=1.0;
const double kHz=1.0e+3*Hz;
const double MHz=1.0e+6*Hz;
const double GHz=1.0e+9*Hz;
const double THz=1.0e+12*Hz;

// Angle
const double radian=1.0;
const double milliradian=1.0e-3*radian;
const double microradian=1.0e-6*radian;
const double degree=0.0174532925199*radian;

// magnetic field
const double Tesla=1.0;
const double kGauss=1.0e-1*Tesla;
const double Gauss=1.0e-3*kGauss;

// mass
const double amu = 931.494028*MeV; // atomic mass unit in MeV

}

const double pi = 4.0*atan(1.0);
const double twoPi = 8.0*atan(1.0);
const double euler = 0.57721566490153286;
