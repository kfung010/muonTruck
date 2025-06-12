#include "detector.hh"

SensitiveDetector::SensitiveDetector(G4String name, G4String detection, EventAction* eventAction) : G4VSensitiveDetector(name), fDetection(detection), fEventAction(eventAction) {}

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
    else if (fDetection == "truckdet") {

        /*FillTruckData(preStepPoint, eventNum, accept);
        if (postStepOnBoundary) {
            FillTruckData(postStepPoint, eventNum, accept);
        }*/
        
        fEventAction->MarkTruckHit(eventNum);               
    }


 
    return true;
}


void SensitiveDetector::FillHitData(G4StepPoint *stepPoint, G4int eventNum) {

    HitData data;
    data.eventNum = eventNum;
    data.muonMomentum = stepPoint->GetMomentum().mag() / GeV;
    data.hitTime = stepPoint->GetGlobalTime();
    data.posMuon = stepPoint->GetPosition();
    data.posDetector = stepPoint->GetTouchable()->GetVolume()->GetTranslation();
    data.momMuon = stepPoint->GetMomentum() / GeV;
    
    fEventAction->CacheRPCHit(eventNum, data);
}


void SensitiveDetector::FillTruckData(G4StepPoint *stepPoint, G4int eventNum) {

    fEventAction->MarkTruckHit(eventNum);
    TruckData data;
    data.eventNum = eventNum;
    data.muonMomentum = stepPoint->GetMomentum().mag() / GeV;
    data.posMuon = stepPoint->GetPosition();
    data.momMuon = stepPoint->GetMomentum() / GeV;
    data.scatTime = stepPoint->GetGlobalTime();
    fEventAction->CacheTruckHit(eventNum, data);        

}








