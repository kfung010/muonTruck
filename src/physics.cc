#include "physics.hh"

PhysicsList::PhysicsList() {
    //RegisterPhysics(new G4EmStandardPhysics_option3());
    //RegisterPhysics(new G4EmStandardPhysics());
    RegisterPhysics(new G4EmStandardPhysics_option4());
    //RegisterPhysics(new G4EmStandardPhysicsWVI());
    //RegisterPhysics(new G4EmStandardPhysicsSS());
    //RegisterPhysics(new G4OpticalPhysics());  
    //RegisterPhysics(new G4DecayPhysics());  
    //RegisterPhysics(new G4RadioactiveDecayPhysics());  
    //RegisterPhysics(new G4StepLimiterPhysics());
}

PhysicsList::~PhysicsList() {}


void PhysicsList::ConstructProcess() {
   
    G4VModularPhysicsList::ConstructProcess();

    G4ParticleDefinition* muon = G4MuonMinus::Definition();
    G4ProcessManager* pmanager = muon->GetProcessManager();
    if (pmanager) {
        pmanager->AddProcess(new G4StepLimiter(), -1, -1, 1); 
    } 
    else {
        G4cerr << "Error: Muon process manager not found!" << G4endl;
    }
}