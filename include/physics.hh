#ifndef PHYSICS_HH
#define PHYSICS_HH
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include <G4EmStandardPhysics_option3.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4EmStandardPhysicsSS.hh>
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysicsWVI.hh"

#include "G4MuonMinus.hh"   
#include "G4ProcessManager.hh" 
#include "G4StepLimiter.hh" 

class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList();
    ~PhysicsList();
    
    virtual void ConstructProcess() override;
};

#endif
