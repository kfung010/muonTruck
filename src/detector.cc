#include "detector.hh"

SensitiveDetector::SensitiveDetector(G4String name, G4String detection) : G4VSensitiveDetector(name), fDetection(detection) {}

SensitiveDetector::~SensitiveDetector(){}

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROIhist) {

    // This function will be triggered whenever a particle enters the physical volume of sensitive detectors
    // If the particle does not have interaction inside the sensitive detectors, then "preStepPoint" records the track information when
    // it just enters the volume, and "postStepPoint" records the track information when it just leaves the volumme
    // But if there are interactions inside the sensitive detectors, a new step is generated for each interaction.
     
    G4Track *track = aStep->GetTrack();

    G4ParticleDefinition* particleDef = track->GetDefinition();
    if (!(particleDef->GetParticleName() == "mu+" || particleDef->GetParticleName() == "mu-")) return true;

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    //G4cout << "preStepPoint height : " << preStepPoint->GetPosition()[1] << G4endl;
    //G4cout << "postStepPoint height : " << postStepPoint->GetPosition()[1] << G4endl;
    
    G4int eventNum = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    if (fDetection == "rpcdet") {
        if (preStepPoint->GetStepStatus() != fGeomBoundary) return true;
        FillHitData(preStepPoint, eventNum);
    }   
    else {
        G4bool preStepOnBoundary = (preStepPoint->GetStepStatus() == fGeomBoundary);
        G4bool postStepOnBoundary = (postStepPoint->GetStepStatus() == fGeomBoundary);
        FillTruckData(preStepPoint, eventNum);
        if (postStepOnBoundary) {
            FillTruckData(postStepPoint, eventNum);
        }
        //FillTruckData(preStepPoint, eventNum);
        //FillTruckData(postStepPoint, eventNum);                
    }

 
    return true;
}


void SensitiveDetector::FillHitData(G4StepPoint *stepPoint, G4int eventNum) {
    G4ThreeVector posMuon = stepPoint->GetPosition();
    G4ThreeVector momMuon = stepPoint->GetMomentum() / GeV;
    const G4VTouchable* touchable = stepPoint->GetTouchable();
    G4VPhysicalVolume* physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
    G4double hitTime = stepPoint->GetGlobalTime();
    G4double muonMomentum = momMuon.mag();

    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager->FillNtupleIColumn(1, 0, eventNum);
    manager->FillNtupleDColumn(1, 1, muonMomentum);
    manager->FillNtupleDColumn(1, 2, hitTime);
    manager->FillNtupleDColumn(1, 3, posMuon.x());
    manager->FillNtupleDColumn(1, 4, posMuon.y());
    manager->FillNtupleDColumn(1, 5, posMuon.z());
    manager->FillNtupleDColumn(1, 6, posDetector.x());
    manager->FillNtupleDColumn(1, 7, posDetector.y());
    manager->FillNtupleDColumn(1, 8, posDetector.z());
    manager->FillNtupleDColumn(1, 9, momMuon.x());
    manager->FillNtupleDColumn(1, 10, momMuon.y());
    manager->FillNtupleDColumn(1, 11, momMuon.z());
    manager->AddNtupleRow(1);
}


void SensitiveDetector::FillTruckData(G4StepPoint *stepPoint, G4int eventNum) {
    G4ThreeVector posMuon = stepPoint->GetPosition();
    G4ThreeVector momMuon = stepPoint->GetMomentum() / GeV;
    G4double hitTime = stepPoint->GetGlobalTime();
    G4double muonMomentum = momMuon.mag();

    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager->FillNtupleIColumn(2, 0, eventNum);
    manager->FillNtupleDColumn(2, 1, muonMomentum);
    manager->FillNtupleDColumn(2, 2, hitTime);
    manager->FillNtupleDColumn(2, 3, posMuon.x());
    manager->FillNtupleDColumn(2, 4, posMuon.y());
    manager->FillNtupleDColumn(2, 5, posMuon.z());
    manager->FillNtupleDColumn(2, 6, momMuon.x());
    manager->FillNtupleDColumn(2, 7, momMuon.y());
    manager->FillNtupleDColumn(2, 8, momMuon.z());
    manager->AddNtupleRow(2);
}













