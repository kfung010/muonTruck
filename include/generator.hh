#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "construction.hh"

class muonGenerator : public G4VUserPrimaryGeneratorAction {
	public:
		muonGenerator();
		~muonGenerator();
		virtual void GeneratePrimaries(G4Event *);

	private:
		G4ParticleGun *fParticleGun;
};

#endif