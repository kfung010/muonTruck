#include "generator.hh"

muonGenerator::muonGenerator(EventAction* eventAction) : fEventAction(eventAction) {

    fMessenger = new G4GenericMessenger(this, "/generation/", "Muon generation");
    fMessenger->DeclareProperty("distribution", distribution, "Muon energy and angular distribution");
    fMessenger->DeclareProperty("muonEnergy", muEne, "Muon energy for monoEnergy_vertical");
    fMessenger->DeclareProperty("muonAngle", muAng, "Muon zenith angle for monoEnergy_tilt");
    fMessenger->DeclareProperty("displacement", muDisp, "Horizontal x displacement from the center");

    G4Random::setTheSeed(12);
    fParticleGun = new G4ParticleGun(1);  //number of particles per event (several events = 1 run)
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName = "mu-";
    G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
    fParticleGun->SetParticleDefinition(particle);
    mass = particle->GetPDGMass();
}

muonGenerator::~muonGenerator() {
    delete fParticleGun;
}

void muonGenerator::GeneratePrimaries(G4Event *anEvent) {

    const RPCConstruction *rpcConstruction = static_cast<const RPCConstruction *> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    G4double xPos, yPos, zPos;
    G4double xMom, yMom, zMom;
    G4double energyGenerate;
    G4double thetaRdmGenerate;

    if (distribution == "monoEnergy_vertical" || distribution == "randomEnergy_vertical") {
        G4double zHighestRPC = rpcConstruction->getzHighestRPC();
        xPos = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
        yPos = G4UniformRand()*rpcConstruction->getyRPCFull() - rpcConstruction->getyRPCFull()/2;
        zPos = zHighestRPC + 10*m;
        xMom = 0;
        yMom = 0;
        zMom = -1;
        if (distribution == "monoEnergy_vertical") energyGenerate = muEne;
        else energyGenerate = G4UniformRand()*9 + 1; 
        thetaRdmGenerate = 0;
    }
    else if (distribution == "monoEnergy_tilt" || distribution == "randomEnergy_tilt") {
        G4double zHighestRPC = rpcConstruction->getzHighestRPC();
        G4double theta = muAng;
        G4double phi = 0;
        xMom = -std::sin(theta)*std::cos(phi);
        yMom = -std::sin(theta)*std::sin(phi);
        zMom = -fabs(std::cos(theta));
        xPos = muDisp*cm-10*m * xMom;
        yPos = -10*m * yMom;
        zPos = zHighestRPC - 10*m * zMom;
        if (distribution == "monoEnergy_tilt") energyGenerate = muEne;
        else energyGenerate = G4UniformRand()*9 + 1; 
        thetaRdmGenerate = 0;
    }
    else if (distribution == "monoEnergy_randomAngle" || distribution == "randomEnergy_randomAngle") {
        G4double xRdmGenerate = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
        G4double yRdmGenerate = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
        G4double zHighestRPC = rpcConstruction->getzHighestRPC();
        G4double thetaRdmGenerate = G4UniformRand()*pi/2;
        G4double phiRdmGenerate = G4UniformRand()*2*pi;
        xMom = -std::sin(thetaRdmGenerate)*std::cos(phiRdmGenerate);
        yMom = -std::sin(thetaRdmGenerate)*std::sin(phiRdmGenerate);
        zMom = -fabs(std::cos(thetaRdmGenerate));
        xPos = xRdmGenerate - 10*m * xMom;
        yPos = yRdmGenerate - 10*m * yMom;
        zPos = zHighestRPC - 10*m * zMom;
        if (distribution == "monoEnergy_randomAngle") energyGenerate = muEne;
        else energyGenerate = G4UniformRand()*49 + 1; 
    }
    else if (distribution == "cosmic") {  
        G4double xRdmGenerate = G4UniformRand()*rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull()/2;
        G4double yRdmGenerate = G4UniformRand()*rpcConstruction->getyRPCFull() - rpcConstruction->getyRPCFull()/2;
        G4double zHighestRPC = rpcConstruction->getzHighestRPC();

        auto thetaPhiEnergy = cosmicMuonZenithPhiAndEnergy();
        thetaRdmGenerate = std::get<0>(thetaPhiEnergy);
        G4double phiRdmGenerate = std::get<1>(thetaPhiEnergy);
        xMom = -std::sin(thetaRdmGenerate)*std::cos(phiRdmGenerate);
        yMom = -std::sin(thetaRdmGenerate)*std::sin(phiRdmGenerate);
        zMom = -fabs(std::cos(thetaRdmGenerate));
        
        xPos = xRdmGenerate - 10*m * xMom;
        yPos = yRdmGenerate - 10*m * yMom;
        zPos = zHighestRPC - 10*m * zMom;
        
        G4double energyRdmGenerate = std::get<2>(thetaPhiEnergy);
        energyGenerate = energyRdmGenerate;
    }
    else throw std::runtime_error("Invalid muon generation mode. Program ended with an error.");

    // Initial position to generate muon
    G4ThreeVector pos(xPos, yPos, zPos);   
    fParticleGun->SetParticlePosition(pos);

    // Initial muon direction
    G4ThreeVector mom(xMom, yMom, zMom); 
    fParticleGun->SetParticleMomentumDirection(mom);

    // Initial muon energy
    //fParticleGun->SetParticleMomentum(momGenerate*GeV);  // This is momentum only. Total energy is p^2+m^2.
    G4double kineticEnergy = energyGenerate * GeV - mass;
    fParticleGun->SetParticleEnergy(kineticEnergy);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    G4int eventNum = anEvent->GetEventID();
    FillGeneratorData(fParticleGun, eventNum);

    //G4ParticleDefinition* particleDef = fParticleGun->GetParticleDefinition();
    //G4ThreeVector position = fParticleGun->GetParticlePosition();
    //G4double momentum = fParticleGun->GetParticleMomentum();
    //G4ThreeVector momentumDirection = fParticleGun->GetParticleMomentumDirection();
    //G4cout << "Particle: " << particleDef->GetParticleName() << G4endl;
    //G4cout << "Momentum: " << momentum << G4endl;
    //G4cout << "Position: (" << position.x() / cm << ", " << position.y() / cm << ", " << position.z() / cm << ") cm" << G4endl;
    //G4cout << "Momentum Direction: (" << momentumDirection.x() << ", " << momentumDirection.y() << ", " << momentumDirection.z() << ")" << G4endl << G4endl;
    
}

void muonGenerator::FillGeneratorData(G4ParticleGun *gun, G4int eventNum) {
    if (fEventAction->GetRecordGenerator()) {
        GeneratorData data;
        data.eventNum = eventNum;
        data.muonEnergy = (gun->GetParticleEnergy() + mass) / GeV;
        data.momMuon = gun->GetParticleMomentumDirection();
        fEventAction->CacheGenerator(eventNum, data);
    }
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
