#include "CollimationParticleGun.h"

/*
CollimationParticleGun::CollimationParticleGun(      
                                const G4String& particleName, 
                                G4double energy,
                                G4ThreeVector position, 
                                G4ThreeVector momentumDirection)
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
*/
CollimationParticleGun::CollimationParticleGun() : G4VUserPrimaryGeneratorAction(), ParticleGun(new G4ParticleGun(1)),  particleTable(G4ParticleTable::GetParticleTable()), particle(nullptr), do_debug(false)
{
	//1 proton at a time
//	ParticleGun  = new G4ParticleGun(1);
//	particleTable = G4ParticleTable::GetParticleTable();
// const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();

/*
	//In the future this could be Pb ions or anything really
	G4String particleName = "proton";

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	particle = particleTable->FindParticle(particleName);

	ParticleGun->SetParticleDefinition(particle);
*/
}

CollimationParticleGun::~CollimationParticleGun()
{
	delete ParticleGun;
}

void CollimationParticleGun::GeneratePrimaries(G4Event* anEvent)
{
	ParticleGun->GeneratePrimaryVertex(anEvent);
}

void CollimationParticleGun::SetParticleDetails(double x, double y, double xp, double yp, double dp, int pdgid, int q)
{
//UNITS MUST BE MeV, mm, rad!
	if(do_debug)
	std::cout << "PDG id: " << pdgid << std::endl;

	if(pdgid < 1000000000)
	{
		particle = particleTable->FindParticle(pdgid);

		if(do_debug)
		std::cout << "name: " << particle->GetParticleName();// << std::endl;
	}
	else
	{
		//Example: 1000822080

		//Now 822080
		int tmpid = pdgid - 1000000000;

		//Now 82208
		tmpid /= 10;

		int A1 = (tmpid%10);
		tmpid /= 10;
		int A2 = (tmpid%10);
		tmpid /= 10;
		int A3 = (tmpid%10);

		int A = A1 + 10*A2 + 100*A3;

		tmpid /= 10;
		int Z1 = (tmpid%10);
		tmpid /= 10;
		int Z2 = (tmpid%10);
		tmpid /= 10;
		int Z3 = (tmpid%10);
		int Z = Z1 + 10*Z2 + 100*Z3;

		if(do_debug)
		std::cout << "Ion: A: " << A << "\tZ: " << Z << "\tQ: " << q; // << std::endl;

		particle = G4IonTable::GetIonTable()->GetIon(pdgid);

	}

	ParticleGun->SetParticleDefinition(particle);

	if(particle->IsGeneralIon())
	{
		ParticleGun->SetParticleCharge(q*CLHEP::eplus);
	}
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

	if(do_debug)
	{
		std::cout << "\tP(GeV): " << p/CLHEP::GeV <<  "\t" << dp/CLHEP::GeV << std::endl;
		std::cout << std::endl;
	}
}

void CollimationParticleGun::SetReferenceEnergy(double energy)
{
	ReferenceEnergy = energy;
}

double CollimationParticleGun::GetReferenceEnergy()
{
	return ReferenceEnergy;
}

void CollimationParticleGun::SetDebug(bool flag)
{
	do_debug = flag;
}

