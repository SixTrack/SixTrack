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
{}

CollimationParticleGun::~CollimationParticleGun()
{
	delete ParticleGun;
}

void CollimationParticleGun::GeneratePrimaries(G4Event* anEvent)
{
	ParticleGun->GeneratePrimaryVertex(anEvent);
}

void CollimationParticleGun::SetParticleDetails(double x, double y, double px, double py, double pz, double e, double p, int pdgid, int q, double mass, double sx, double sy, double sz)
{
//UNITS MUST BE MeV, mm, rad!
	if(do_debug)
	std::cout << "GEANT4> PDG id: " << pdgid << std::endl;

	if(pdgid < 1000000000)
	{
		particle = particleTable->FindParticle(pdgid);

		if(do_debug)
		std::cout << "GEANT4> name: " << particle->GetParticleName();// << std::endl;
	}
	else
	{
		//Example: 1000822080

		if(do_debug)
		{
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

			std::cout << "GEANT4> Ion: A: " << A << "\tZ: " << Z << "\tQ: " << q; // << std::endl;
		}

		//Remove states other than the ground state otherwise we get a crash.
		while(pdgid%10 !=0)
		{
			pdgid--;
		}

		particle = G4IonTable::GetIonTable()->GetIon(pdgid);

	}

	ParticleGun->SetParticleDefinition(particle);

	if(particle->IsGeneralIon())
	{
		ParticleGun->SetParticleCharge(q*CLHEP::eplus);
	}
	const G4double mp = particle->GetPDGMass();
	if(mp != mass)
	{
		std::cout.precision(17);
		std::cout << "GEANT4> Mass missmatch between Geant4 and SixTrack!" << std::endl;
		std::cout << "GEANT4> PDG mass (MeV)     : " << mp << std::endl;
		std::cout << "GEANT4> SixTrack mass (MeV): " << mass << std::endl;
		std::cout << "GEANT4> Please set the particle mass in SixTrack to the Geant4 (PDG) value!" << std::endl;
		std::cout << "GEANT4> If the HION block is used, convert to GeV!" << std::endl;
		exit(EXIT_FAILURE);
	}

	//The kinetic energy (Total energy - rest mass)
	ParticleGun->SetParticleEnergy(e - mp);

	//How to deal with longitudinal coordinate?
	ParticleGun->SetParticlePosition(G4ThreeVector(x, y, 0));
	ParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
//	ParticleGun->SetParticleMomentum(G4ThreeVector(px,py,pz));

	ParticleGun->SetParticlePolarization(G4ThreeVector(sx,sy,sz));

	if(do_debug)
	{
		std::cout << "\tP(GeV): " << p/CLHEP::GeV <<  "\t" << e/CLHEP::GeV << std::endl;
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

