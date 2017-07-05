#include "CollimationParticleGun.h"

/*
CollimationParticleGun::CollimationParticleGun(      
                                const G4String& particleName, 
                                G4double energy,
                                G4ThreeVector position, 
                                G4ThreeVector momentumDirection)
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
*/
CollimationParticleGun::CollimationParticleGun() : G4VUserPrimaryGeneratorAction(), ParticleGun(nullptr), particle(nullptr)
{
	//1 proton at a time
	ParticleGun  = new G4ParticleGun(1);

	//In the future this could be Pb ions or anything really
	G4String particleName = "proton";

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	particle = particleTable->FindParticle(particleName);

	ParticleGun->SetParticleDefinition(particle);
}

CollimationParticleGun::~CollimationParticleGun()
{
	delete ParticleGun;
}

void CollimationParticleGun::GeneratePrimaries(G4Event* anEvent)
{
	ParticleGun->GeneratePrimaryVertex(anEvent);
}

void CollimationParticleGun::SetParticleDetails(double x, double y, double xp, double yp, double dp)
{
//UNITS MUST BE MeV, mm, rad!

	const G4double mp = particle->GetPDGMass();
/*
	G4double BeamMomentum = sqrt((ReferenceEnergy*ReferenceEnergy) - (mp*mp));
	G4double ParticleMomentum = BeamMomentum * (1.0 + dp);
	G4double ParticleKE = sqrt((ParticleMomentum*ParticleMomentum) + (mp*mp)) - mp;
	G4cout << ParticleKE << "\t" << ParticleMomentum << "\t" << BeamMomentum << "\t" << ReferenceEnergy << "\t" << dp << std::endl;
*/

	//The kinetic energy (Total energy - rest mass)
	ParticleGun->SetParticleEnergy(dp - mp);
	double p = sqrt(dp*dp - mp*mp);
	//How to deal with longitudinal coordinate?
	ParticleGun->SetParticlePosition(G4ThreeVector(x, y, 0));
	ParticleGun->SetParticleMomentumDirection(G4ThreeVector(xp,yp,p));
}

void CollimationParticleGun::SetReferenceEnergy(double energy)
{
	ReferenceEnergy = energy;
}

double CollimationParticleGun::GetReferenceEnergy()
{
	return ReferenceEnergy;
}
