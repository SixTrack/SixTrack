/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
//
/////////////////////////////////////////////////////////////////////////

#ifndef RandomNG_h
#define RandomNG_h 1

//#include "merlin_config.h"
#include <cassert>

class ACG;
class Normal;
class Uniform;
//class Poisson;
class Landau;

/**
* A class which represents a single random number
* generator.
*/
class RandGenerator
{
public:

	RandGenerator (unsigned  iseed = 0);
	~RandGenerator ();

	unsigned  getSeed () const;

	/**
	* Resets the seed for the generators to the last supplied
	* seed value.
	*/
	void reset ();

	/**
	* Resets the initial seed for the generator.
	*/
	void reset (unsigned  iseed);

	/**
	* Generates a random number from a uniform (Gaussian)
	* distribution with the specified mean and variance.
	*/
	double normal (double mean, double variance);

	/**
	* Generates a random number from a uniform (Gaussian)
	* distribution with the specified mean and variance. The
	* resulting distribution is truncated to +/-cutoff
	* standard deviations.
	*/
	double normal (double mean, double variance, double cutoff);

	/**
	* Generates a uniform random number in the range
	* |low,high> inclusive.
	*/
	double uniform (double low, double high);

	/**
	* Generates a poisson random number.
	*/
	//double poisson (double u);

	/**
	* Generates a landau random number.
	*/
	double landau ();

	/**
	* Initialised the generator. Should be called before any
	* other generator function.
	*/
	void init (unsigned  iseed = 0);

private:

	unsigned  nseed;

	ACG* gen;
	Normal* gaussGen;
	Uniform* uniformGen;
	//Poisson* poissonGen;
	Landau* landauGen;

	void ResetGenerators ();

	RandGenerator(const RandGenerator& rand);
	RandGenerator& operator=(const RandGenerator& rand);

};

/**
* Singleton class for generating continuous floating point
* numbers from specific distributions. Currently only the
* normal (Gaussian) and uniform distributions are
* supported.
*/
class RandomNG
{
public:
	/*
	* Resets the seed for the generators to the last supplied
	* seed value.
	*/
	static void reset ();

	/**
	* Resets the initial seed for the generator.
	*/
	static void reset (unsigned  iseed);

	/**
	* Get the seed for the generator.
	*/
	static unsigned getSeed();

	/**
	* Generates a random number from a uniform (Gaussian)
	* distribution with the specified mean and variance.
	*/
	static double normal (double mean, double variance);

	/**
	* Generates a random number from a uniform (Gaussian)
	* distribution with the specified mean and variance. The
	* resulting distribution is truncated to +/-cutoff
	* standard deviations.
	*/
	static double normal (double mean, double variance, double cutoff);

	/**
	* Generates a uniform random number in the range
	* |low,high> inclusive.
	*/
	static double uniform (double low, double high);

	/**
	* Generates a uniform random number in the range
	* |low,high> inclusive.
	*/
	//static double poisson (double u);

	/**
	* Generates a uniform random number in the range
	* |low,high> inclusive.
	*/
	static double landau ();

	/**
	* Initialised the generator. Should be called before any
	* other generator function.
	*/
	static void init (unsigned  iseed = 0);

private:

	static RandGenerator* generator;
};

inline unsigned RandGenerator::getSeed () const
{
	return nseed;
}

inline unsigned RandomNG::getSeed ()
{
	assert(generator);
	return generator->getSeed();
}
inline void RandomNG::reset ()
{
	assert(generator);
	generator->reset();
}

inline void RandomNG::reset (unsigned iseed)
{
	assert(generator);
	generator->reset(iseed);
}

inline double RandomNG::normal (double mean, double variance)
{
	assert(generator);
	return generator->normal(mean,variance);
}

inline double RandomNG::normal (double mean, double variance, double cutoff)
{
	assert(generator);
	return generator->normal(mean,variance,cutoff);
}

inline double RandomNG::uniform (double low, double high)
{
	assert(generator);
	return generator->uniform(low,high);
}
/*
inline double RandomNG::poisson (double u)
{
	assert(generator);
	return generator->poisson(u);
}
*/
inline double RandomNG::landau ()
{
	assert(generator);
	return generator->landau();
}

inline void RandomNG::init (unsigned iseed)
{
	if(generator)
	{
		delete generator;
	}
	generator = new RandGenerator(iseed);
}

#endif
