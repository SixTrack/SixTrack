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
	void SetParticleDetails(double x, double y, double xp, double yp, double dp, int pdgid, int q);
	void SetReferenceEnergy(double);
	double GetReferenceEnergy();

private:
	G4ParticleGun* ParticleGun;
	G4ParticleTable* particleTable;
//	const G4IonTable* ionTable;
	G4ParticleDefinition* particle;
	G4double ReferenceEnergy;
};

#endif

