#ifndef _Landau_h_
#define _Landau_h_
#include "Random.h"

class Landau: public Random
{
//const float F[982];

public:
	Landau(RNG *gen);

	virtual double operator()();
};


inline Landau::Landau(RNG *gen) : Random(gen) {}

#endif
