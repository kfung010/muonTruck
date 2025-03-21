#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "G4PhysicalConstants.hh"
#include "construction.hh"

class muonGenerator : public G4VUserPrimaryGeneratorAction {
    public:
        muonGenerator();
        ~muonGenerator();
        virtual void GeneratePrimaries(G4Event *);

    private:
        G4ParticleGun *fParticleGun;
        G4double cosThetaStar(G4double theta);
        G4double flux(G4double theta, G4double ene);
        std::tuple<G4double, G4double, G4double> cosmicMuonZenithPhiAndEnergy(G4double thetaLow=0, G4double thetaUp=pi/2, G4double phiLow=0, G4double phiUp=2*pi, G4double eneLow=0.5, G4double eneUp=10);
		
        G4String distribution;
        G4GenericMessenger *fMessenger;
};

#endif