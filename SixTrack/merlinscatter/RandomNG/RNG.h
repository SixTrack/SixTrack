// This may look like C code, but it is really -*- C++ -*-
/*
Copyright (C) 1988 Free Software Foundation
    written by Dirk Grunwald (grunwald@cs.uiuc.edu)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#ifndef _RNG_h
#define _RNG_h 1
#ifdef __GNUG__
#endif

#include <cassert>
#include <cmath>
//#include <_G_config.h>

// 32 bit integer types
typedef unsigned int _G_uint32_t;
typedef int _G_int32_t;

union PrivateRNGSingleType  		   	// used to access floats as unsigneds
{
	float s;
	unsigned int u;
};

union PrivateRNGDoubleType  		   	// used to access doubles as unsigneds
{
	double d;
	unsigned int u[2];
};

//
// Base class for Random Number Generators. See ACG and MLCG for instances.
//
class RNG
{
	static PrivateRNGSingleType singleMantissa;	// mantissa bit vector
	static PrivateRNGDoubleType doubleMantissa;	// mantissa bit vector
public:
	RNG();
	virtual ~RNG() {};
	//
	// Return a long-words word of random bits
	//
	virtual _G_uint32_t asLong() = 0;
	virtual void reset() = 0;
	//
	// Return random bits converted to either a float or a double
	//
	float asFloat();
	double asDouble();
};

#endif
