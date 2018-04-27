#ifndef _helperfunctions_h_
#define _helperfunctions_h_

#include <string>
#include <vector>
#include <complex>

#include "RandomNG/RandomNG.h"

typedef std::complex<double> Complex;

/**
* Rather than have thousands of files moved over from Merlin, we are just going to copy all the required classes into these files
*/

/**
* Root class for all Merlin exceptions.
*/
class MerlinException
{
public:

	/**
	* Constructor: Builds the MerlinException and sets the exception message.
	* @param[in] s The exception message.
	*/
	explicit MerlinException(const std::string& s);

	/**
	* Constructor: Builds the MerlinException.
	*/
	MerlinException();

	/**
	* Virtual destructor.
	*/
	virtual ~MerlinException() {};

	/**
	* Gets the exception message.
	* @return The current exception message.
	*/
	const char* Msg() const;

protected:

	/**
	* Sets the exception message.
	* @see AppendMsg
	* @param[in] s The exception message.
	*/
	void SetMsg(const std::string& s);

	/**
	* Appends a string to the exception message.
	* @see SetMsg
	* @param[in] s The string to append to the exception message.
	*/
	void AppendMsg(const std::string& s);

private:

	/**
	* String storage containing the exception message.
	*/
	std::string msg;
};

inline MerlinException::MerlinException(const std::string& s)
	: msg(s)
{}

inline MerlinException::MerlinException()
	: msg()
{}

inline const char* MerlinException::Msg() const
{
	return msg.c_str();
}

inline void MerlinException::SetMsg(const std::string& s)
{
	msg=s;
}

inline void MerlinException::AppendMsg(const std::string& s)
{
	msg+=s;
}

class RangeBase
{
public:
	typedef enum {ok,belowlower,aboveupper} Result;
	virtual ~RangeBase() {};
};

//	Represents an allowed contiguous range of a floating
//	point number.
struct UnboundedRange {};

template <class T, class C = std::less<T> >
class NumericalRange : public RangeBase
{
public:

	//	Constructor taking the range (lo,hi>
	NumericalRange (const T& lo, const T& hi);

	//    Construct an unbounded range
	NumericalRange (UnboundedRange) : lower(1),upper(-1) {}

	//    Construct a fixed point (zero range)
	explicit NumericalRange (const T& fp) : lower(fp),upper(fp) {}

	//	Returns true if this range is unbouned (i.e. represents
	//	+/-infinity).
	bool IsUnbounded () const;
	bool IsFixedPoint () const
	{
		return lower==upper;
	}

	//	Returns true if lower<=x<=upper.
	bool operator () (const T& x) const;
	bool operator == (const NumericalRange& rhs) const;
	bool operator != (const NumericalRange& rhs) const;

	//	Checks id x is within the range. Returns belowlower, ok or
	//	aboveupper.
	RangeBase::Result Check (const T& x) const;

	//	The lowerimum value of the range.
	T lower;

	//	The upperimum value of the range.
	T upper;
};


template <class T, class C>
inline NumericalRange<T,C>::NumericalRange (const T& lo, const T& hi)
	: lower(lo),upper(hi)
{
	if(C()(upper,lower))
	{
		std::swap(lower,upper);
	}
}

template <class T, class C>
inline bool NumericalRange<T,C>::IsUnbounded () const
{
	return fequal(lower,1.0) && fequal(upper,-1.0);
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator () (const T& x) const
{
	if(IsUnbounded())
	{
		return true;
	}
	else
	{
		return x>=lower && x<=upper;
	}
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator == (const NumericalRange& rhs) const
{
	return lower==rhs.lower && upper==rhs.upper;
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator != (const NumericalRange& rhs) const
{
	return lower!=rhs.lower || upper!=rhs.upper;
}

template <class T, class C>
inline RangeBase::Result NumericalRange<T,C>::Check (const T& x) const
{
//dk    if(isUnbounded())
	if(IsUnbounded())
	{
		return ok;
	}
	else if(x<lower)
	{
		return belowlower;
	}
	else if(x<=upper)
	{
		return ok;
	}
	else
	{
		return aboveupper;
	}
}

// float numerical range
typedef NumericalRange< double > FloatRange;

/**
* class Interpolation
* An interpolation functor which interpolates values from a data table.
* Currenly only linear interpolation is assumed.
*/
class Interpolation
{
public:

	// exception
	class BadRange : public MerlinException
	{
	public:
		BadRange(double x, const FloatRange& r);
		double value;
		FloatRange valid_range;
	};

	// Implementation method for interpolation
	class Method
	{
	public:
		virtual ~Method() {}
		virtual double ValueAt(double x) const =0;
	};

	// Interpolation of equally spaced data points
	Interpolation(const std::vector<double>& yvals, double xmin, double dx);

	// Interpolation of arbitrary spaced data points
	Interpolation(const std::vector<double>& xvals, const std::vector<double>& yvals);

	~Interpolation();

	double operator()(double x) const
	{
		return itsMethod->ValueAt(x);
	}

private:

	Method* itsMethod;
};


namespace PhysicalConstants
{
// Data Members for Class Attributes

extern const double Avogadro;

extern const double AtomicMassUnit;

extern const double SpeedOfLight;

extern const double ElectronMass;

extern const double ProtonMass;

extern const double PlanckConstant;

extern const double PlanckConstantBar;

extern const double ElectronCharge;

extern const double ElectronMassMeV;

extern const double ProtonMassMeV;

extern const double ProtonMassGeV;

extern const double FreeSpacePermeability;

extern const double FreeSpacePermittivity;

extern const double ElectronRadius;

extern const double ElectronGe;

extern const double MuonMass;

extern const double MuonMassMeV;

extern const double MuonLifetime;

extern const double FineStructureConstant;

extern const double PionZeroMassMeV;

extern double LorentzBeta (double gamma);

extern double LorentzGamma (double beta);
extern double LorentzGamma (double momentum, double mass);
}

// Merlin common units

namespace PhysicalUnits
{

//Area
extern const double barn;

// length
extern const double meter;
extern const double centimeter;
extern const double millimeter;
extern const double micrometer;
extern const double nanometer;

// time
extern const double second;
extern const double millisecond;
extern const double microsecond;
extern const double nanosecond;
extern const double picosecond;
extern const double minute;
extern const double hour;
extern const double day;
extern const double year;

// energy
extern const double eV;
extern const double keV;
extern const double MeV;
extern const double GeV;
extern const double TeV;

// voltage
extern const double Volt;
extern const double kV;
extern const double MV;
extern const double GV;
extern const double TV;

// Frequency
extern const double Hz;
extern const double kHz;
extern const double MHz;
extern const double GHz;
extern const double THz;

// Angle
extern const double radian;
extern const double milliradian;
extern const double microradian;
extern const double degree;

// magnetic field
extern const double Tesla;
extern const double kGauss;
extern const double Gauss;

// mass
extern const double amu;

}

extern const double pi;
extern const double twoPi;
extern const double euler;

#endif
