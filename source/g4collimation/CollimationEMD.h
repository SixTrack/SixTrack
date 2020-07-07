//Provided by L.Nevay (RHUL) via BDSIM
//Needed to correctly simulate heavy ion collimation
#ifndef CollimationEMD_h_ 
#define CollimationEMD_h_

#include "G4VPhysicsConstructor.hh"

class EMDissociation : public G4VPhysicsConstructor
{
public:
	EMDissociation();
	virtual void ConstructParticle();
	virtual void ConstructProcess();
};

#endif

