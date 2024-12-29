#include "generator.hh"

muonGenerator::muonGenerator() {
    fParticleGun = new G4ParticleGun(1);  //number of particles per event (several events = 1 run)
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "mu-";
    G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);
}
muonGenerator::~muonGenerator() {
    delete fParticleGun;
}

void muonGenerator::GeneratePrimaries(G4Event *anEvent) {

	const RPCConstruction *rpcConstruction = static_cast<const RPCConstruction *> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	const double PI = 3.14159265358979323846;
	
	G4double xRdmGenerate = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
	G4double yHighestRPC = rpcConstruction->getyHighestRPC();
	G4double zRdmGenerate = G4UniformRand()*rpcConstruction->getzRPCFull() - rpcConstruction->getzRPCFull()/2;
	
	G4double thetaRdmGenerate = G4UniformRand()*PI/4;
	G4double phiRdmGenerate = G4UniformRand()*2*PI;
	G4double xMomRdmGenerate = std::sin(thetaRdmGenerate)*std::cos(phiRdmGenerate);
	G4double yMomRdmGenerate = fabs(std::cos(thetaRdmGenerate));
	G4double zMomRdmGenerate = std::sin(thetaRdmGenerate)*std::sin(phiRdmGenerate);

	G4double xPos = xRdmGenerate + 10*m * xMomRdmGenerate;
	G4double yPos = yHighestRPC + 10*m * yMomRdmGenerate;
	G4double zPos = zRdmGenerate + 10*m * zMomRdmGenerate;

	G4ThreeVector pos(xPos, yPos, zPos);   
	fParticleGun->SetParticlePosition(pos);

	G4ThreeVector mom(-xMomRdmGenerate, -yMomRdmGenerate, -zMomRdmGenerate); 
    fParticleGun->SetParticleMomentumDirection(mom);

	G4double momGenerate = G4UniformRand()*9 + 1;
    fParticleGun->SetParticleMomentum(momGenerate*GeV);  // This is momentum only. Total energy is p^2+m^2.

    fParticleGun->GeneratePrimaryVertex(anEvent);

	G4ParticleDefinition* particleDef = fParticleGun->GetParticleDefinition();

	G4ThreeVector position = fParticleGun->GetParticlePosition();
	G4double momentum = fParticleGun->GetParticleMomentum();
	G4ThreeVector momentumDirection = fParticleGun->GetParticleMomentumDirection();
	//G4cout << "Particle: " << particleDef->GetParticleName() << G4endl;
	//G4cout << "Momentum: " << momentum << G4endl;
	//G4cout << "Position: (" << position.x() / cm << ", " << position.y() / cm << ", " << position.z() / cm << ") cm" << G4endl;
	//G4cout << "Momentum Direction: (" << momentumDirection.x() << ", " << momentumDirection.y() << ", " << momentumDirection.z() << ")" << G4endl << G4endl;

	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	manager->FillNtupleDColumn(0, 0, momentum);
	manager->AddNtupleRow(0);
}