//Provided by L.Nevay (RHUL) via BDSIM
//Needed to correctly simulate heavy ion collimation

#include "CollimationEMD.h"

#include "G4EMDissociation.hh"
#include "G4EMDissociationCrossSection.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4IonConstructor.hh"
#include "G4ProcessManager.hh"

EMDissociation::EMDissociation() : G4VPhysicsConstructor("EMDissociation") {}

void EMDissociation::ConstructParticle()
{
  G4Gamma::Gamma();
  G4GenericIon::GenericIon();

  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void EMDissociation::ConstructProcess()
{
  G4HadronInelasticProcess* inelProcIon = new G4HadronInelasticProcess("ionInelastic", G4GenericIon::GenericIon());

  G4EMDissociationCrossSection* crossSectionData = new G4EMDissociationCrossSection();
  inelProcIon->AddDataSet(crossSectionData);

  G4EMDissociation* emdModel = new G4EMDissociation();
  emdModel->SetMaxEnergy(100*CLHEP::TeV);
  inelProcIon->RegisterMe(emdModel);

  G4ProcessManager* pmanager = G4GenericIon::GenericIon()->GetProcessManager();
  pmanager->AddDiscreteProcess(inelProcIon);
  
}

