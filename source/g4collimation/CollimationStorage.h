#ifndef _g4_collimation_storage_
#define _g4_collimation_storage_

#include "G4Types.hh"

//g4_add_particle(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nzz(j), naa(j), nqq(j))
struct G4Stuff
{
	G4double x;
	G4double y;
	G4double px;
	G4double py;
	G4double e;
	G4double p;
	G4double m;
	G4int pdgid;
	G4int a;
	G4int z;
	G4int q;
	G4int id;
};

#endif

