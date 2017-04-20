// This may look like C code, but it is really -*- C++ -*-
/*
Copyright (C) 1989 Free Software Foundation

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
#ifdef __GNUG__
#endif
#include "MLCG.h"
//
//	SEED_TABLE_SIZE must be a power of 2
//


#define SEED_TABLE_SIZE 32

static _G_int32_t seedTable[SEED_TABLE_SIZE] =
{
	_G_int32_t(0xbdcc47e5), _G_int32_t(0x54aea45d), _G_int32_t(0xec0df859), _G_int32_t(0xda84637b),
	_G_int32_t(0xc8c6cb4f), _G_int32_t(0x35574b01), _G_int32_t(0x28260b7d), _G_int32_t(0x0d07fdbf),
	_G_int32_t(0x9faaeeb0), _G_int32_t(0x613dd169), _G_int32_t(0x5ce2d818), _G_int32_t(0x85b9e706),
	_G_int32_t(0xab2469db), _G_int32_t(0xda02b0dc), _G_int32_t(0x45c60d6e), _G_int32_t(0xffe49d10),
	_G_int32_t(0x7224fea3), _G_int32_t(0xf9684fc9), _G_int32_t(0xfc7ee074), _G_int32_t(0x326ce92a),
	_G_int32_t(0x366d13b5), _G_int32_t(0x17aaa731), _G_int32_t(0xeb83a675), _G_int32_t(0x7781cb32),
	_G_int32_t(0x4ec7c92d), _G_int32_t(0x7f187521), _G_int32_t(0x2cf346b4), _G_int32_t(0xad13310f),
	_G_int32_t(0xb89cff2b), _G_int32_t(0x12164de1), _G_int32_t(0xa865168d), _G_int32_t(0x32b56cdf)
};

MLCG::MLCG(_G_int32_t seed1, _G_int32_t seed2)
{
	initialSeedOne = seed1;
	initialSeedTwo = seed2;
	reset();
}

void
MLCG::reset()
{
	_G_int32_t seed1 = initialSeedOne;
	_G_int32_t seed2 = initialSeedTwo;

	//
	//	Most people pick stupid seed numbers that do not have enough
	//	bits. In this case, if they pick a small seed number, we
	//	map that to a specific seed.
	//
	if (seed1 < 0)
	{
		seed1 = (seed1 + 2147483561);
		seed1 = (seed1 < 0) ? -seed1 : seed1;
	}

	if (seed2 < 0)
	{
		seed2 = (seed2 + 2147483561);
		seed2 = (seed2 < 0) ? -seed2 : seed2;
	}

	if (seed1 > -1 && seed1 < SEED_TABLE_SIZE)
	{
		seedOne = seedTable[seed1];
	}
	else
	{
		seedOne = seed1 ^ seedTable[seed1 & (SEED_TABLE_SIZE-1)];
	}

	if (seed2 > -1 && seed2 < SEED_TABLE_SIZE)
	{
		seedTwo = seedTable[seed2];
	}
	else
	{
		seedTwo = seed2 ^ seedTable[ seed2 & (SEED_TABLE_SIZE-1) ];
	}
	seedOne = (seedOne % 2147483561) + 1;
	seedTwo = (seedTwo % 2147483397) + 1;
}

_G_uint32_t MLCG::asLong()
{
	_G_int32_t k = seedOne % 53668;

	seedOne = 40014 * (seedOne-k * 53668) - k * 12211;
	if (seedOne < 0)
	{
		seedOne += 2147483563;
	}

	k = seedTwo % 52774;
	seedTwo = 40692 * (seedTwo - k * 52774) - k * 3791;
	if (seedTwo < 0)
	{
		seedTwo += 2147483399;
	}

	_G_int32_t z = seedOne - seedTwo;
	if (z < 1)
	{
		z += 2147483562;
	}
	return( static_cast<unsigned long>(z) );
}

