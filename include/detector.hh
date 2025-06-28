#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "eventAction.hh"

class SensitiveDetector: public G4VSensitiveDetector {

    public:
        SensitiveDetector(G4String name, EventAction* eventAction);
        ~SensitiveDetector();
        
        enum DetectionType {
            RPC_DETECTION,
            BOX_DETECTION,
            WORLD_DETECTION
        };
        
    private:
        DetectionType fDetectionType;
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
        void FillHitData(G4StepPoint* stepPoint, G4int eventNum);
        void FillScatteringData(G4StepPoint* stepPoint, G4int eventNum);
        EventAction *fEventAction;
};


#endif