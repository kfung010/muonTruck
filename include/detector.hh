#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsOrderedFreeVector.hh"

class SensitiveDetector: public G4VSensitiveDetector {

    public:
        SensitiveDetector(G4String name, G4String detection);
        ~SensitiveDetector();
    private:
        G4String fDetection;
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
        void FillHitData(G4StepPoint* stepPoint, G4int eventNum);
        void FillTruckData(G4StepPoint* stepPoint, G4int eventNum);
};

#endif