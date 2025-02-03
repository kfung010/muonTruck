#include "generator.hh"

muonGenerator::muonGenerator() {
	G4Random::setTheSeed(12);
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
	
	// Generate muon vertically above the center
	/*G4double yHighestRPC = rpcConstruction->getyHighestRPC();
	G4double xPos = 0 + 10*m * 0;
	G4double yPos = yHighestRPC + 10*m * 1;
	G4double zPos = 0 + 10*m * 0;
	G4double xMom = 0;
	G4double yMom = -1;
	G4double zMom = 0;
	G4double momGenerate = 0.5;*/

	// Generate muon from 0 to 45 degrees to the vertical, and it is ensured to pass through the top RPC plate
	/*G4double xRdmGenerate = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
	G4double yHighestRPC = rpcConstruction->getyHighestRPC();
	G4double zRdmGenerate = G4UniformRand()*rpcConstruction->getzRPCFull() - rpcConstruction->getzRPCFull()/2;
	G4double thetaRdmGenerate = G4UniformRand()*pi/4;
	G4double phiRdmGenerate = G4UniformRand()*2*pi;
	G4double xMom = -std::sin(thetaRdmGenerate)*std::cos(phiRdmGenerate);
	G4double yMom = -fabs(std::cos(thetaRdmGenerate));
	G4double zMom = -std::sin(thetaRdmGenerate)*std::sin(phiRdmGenerate);
	G4double xPos = xRdmGenerate - 10*m * xMom;
	G4double yPos = yHighestRPC - 10*m * yMom;
	G4double zPos = zRdmGenerate - 10*m * zMom;
	G4double momGenerate = G4UniformRand()*9 + 1;*/

	// Generate muon energy and direction according to a measured distribution, and it is ensured to pass through the center of the top RPC plate
	G4double xCenter = 0;
	G4double yHighestRPC = rpcConstruction->getyHighestRPC();
	G4double zCenter = 0;
	auto thetaPhiEnergy = cosmicMuonZenithPhiAndEnergy();
	G4double thetaRdmGenerate = std::get<0>(thetaPhiEnergy);
	G4double phiRdmGenerate = std::get<1>(thetaPhiEnergy);
	G4double xMom = -std::sin(thetaRdmGenerate)*std::cos(phiRdmGenerate);
	G4double yMom = -fabs(std::cos(thetaRdmGenerate));
	G4double zMom = -std::sin(thetaRdmGenerate)*std::sin(phiRdmGenerate);
	G4double xPos = xCenter - 10*m * xMom;
	G4double yPos = yHighestRPC - 10*m * yMom;
	G4double zPos = zCenter - 10*m * zMom;
	G4double energyRdmGenerate = std::get<2>(thetaPhiEnergy);
	G4double momGenerate = sqrt(energyRdmGenerate*energyRdmGenerate-0.105*0.105);

	// Initial position to generate muon
	G4ThreeVector pos(xPos, yPos, zPos);   
	fParticleGun->SetParticlePosition(pos);

	// Initial muon direction
	G4ThreeVector mom(xMom, yMom, zMom); 
    fParticleGun->SetParticleMomentumDirection(mom);

	// Initial muon energy
    fParticleGun->SetParticleMomentum(momGenerate*GeV);  // This is momentum only. Total energy is p^2+m^2.


    fParticleGun->GeneratePrimaryVertex(anEvent);

	G4ParticleDefinition* particleDef = fParticleGun->GetParticleDefinition();

	// Store the generation information in the ntuple
	G4ThreeVector position = fParticleGun->GetParticlePosition();
	G4double momentum = fParticleGun->GetParticleMomentum();
	G4ThreeVector momentumDirection = fParticleGun->GetParticleMomentumDirection();
	//G4cout << "Particle: " << particleDef->GetParticleName() << G4endl;
	//G4cout << "Momentum: " << momentum << G4endl;
	//G4cout << "Position: (" << position.x() / cm << ", " << position.y() / cm << ", " << position.z() / cm << ") cm" << G4endl;
	//G4cout << "Momentum Direction: (" << momentumDirection.x() << ", " << momentumDirection.y() << ", " << momentumDirection.z() << ")" << G4endl << G4endl;

	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	manager->FillNtupleDColumn(0, 0, momentum);
	manager->FillNtupleDColumn(0, 1, momGenerate*xMom);
	manager->FillNtupleDColumn(0, 2, momGenerate*yMom);
	manager->FillNtupleDColumn(0, 3, momGenerate*zMom);
	manager->FillNtupleDColumn(0, 4, thetaRdmGenerate);
	manager->AddNtupleRow(0);
}

G4double muonGenerator::cosThetaStar(G4double theta) {
	G4double P1 = 0.102573;
	G4double P2 = -0.068287;
	G4double P3 = 0.958633;
	G4double P4 = 0.0407253;
	G4double P5 = 0.817285;
	G4double numer = cos(theta)*cos(theta) + P1*P1 + P2*pow(cos(theta), P3) + P4*pow(cos(theta), P5);
	G4double denom = 1 + P1*P1 + P2 + P4;
	return sqrt(numer/denom);
}

G4double muonGenerator::flux(G4double theta, G4double ene) {
    G4double part1 = 0.14;
	G4double part2 = pow(ene*(1+3.64/ene/pow(cosThetaStar(theta), 1.29)),-2.7);
	G4double part3 = 1/(1+1.1*ene*cosThetaStar(theta)/115) + 0.054/(1+1.1*ene*cosThetaStar(theta)/850);
	return part1*part2*part3;
}

std::tuple<G4double, G4double, G4double> muonGenerator::cosmicMuonZenithPhiAndEnergy(G4double thetaLow, G4double thetaUp, G4double phiLow, G4double phiUp, G4double eneLow, G4double eneUp) {

	G4double theta, ene, phi, u;

	G4double f_max = 10;  //Equal to or larger than the max value of flux()
	while (true) {
		theta = G4UniformRand()*(thetaUp-thetaLow) + thetaLow;
		phi = G4UniformRand()*(phiUp-phiLow) + phiLow;
		ene = G4UniformRand()*(eneUp-eneLow) + eneLow;
		u = G4UniformRand();
        if (u <= flux(theta, ene) / f_max)
			break;
    }

	return std::make_tuple(theta, phi, ene);
}


