#ifndef CollimationParticleGun_h
#define CollimationParticleGun_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class CollimationParticleGun: public G4VUserPrimaryGeneratorAction
{
public:
	CollimationParticleGun();
	~CollimationParticleGun();

	void GeneratePrimaries(G4Event*);
	void SetParticleDetails(double x, double y, double xp, double yp, double zp, double e, double p, int pdgid, int q, double mass, double sx, double sy, double sz);
	void SetReferenceEnergy(double);
	double GetReferenceEnergy();
	void SetDebug(bool);

private:
	G4ParticleGun* ParticleGun;
	G4ParticleTable* particleTable;
	G4ParticleDefinition* particle;
	G4double ReferenceEnergy;
	bool do_debug;
};

#endif

