#ifndef _g4_collimation_storage_
#define _g4_collimation_storage_

#include <string>

#include "G4Types.hh"

//g4_add_particle(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nzz(j), naa(j), nqq(j))
struct CollimationParticle
{
	G4double x;
	G4double y;
	G4double px;
	G4double py;
	G4double pz;
	G4double e;
	G4double p;
	G4double m;
	G4double sx;
	G4double sy;
	G4double sz;
	G4double t;
	G4int pdgid;
	G4int a;
	G4int z;
	G4int q;
	G4int id;
	G4int parent_id;
};

struct CollimationEnergyDeposition
{
	std::string name;
	G4int type;
	G4double xstep;
	G4double ystep;
	G4double zstep;

	G4double xmax;
	G4double ymax;
	G4double zmax;

	G4double xbigstep;
	G4double ybigstep;
	G4double zbigstep;

	G4double xbigmax;
	G4double ybigmax;
	G4double zbigmax;

};

struct find_edep_name
{
	std::string fname;
	find_edep_name(std::string n) : fname(n) {}
	bool operator()(const CollimationEnergyDeposition& l)
	{
		return l.name == fname;
	}
};

//#include <algorithm>
//std::vector<CollimationEnergyDeposition>::iterator itr = std::find_if(var.begin(), var.end(), find_edep_name(name));

#endif

