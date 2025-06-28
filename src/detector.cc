#include "detector.hh"

SensitiveDetector::SensitiveDetector(G4String name, EventAction* eventAction) : G4VSensitiveDetector(name), fEventAction(eventAction) {

    if (name.find("rpcdet") != G4String::npos) {
        fDetectionType = RPC_DETECTION;
    } 
    else if (name.find("worlddet") != G4String::npos) {
        fDetectionType = WORLD_DETECTION;
    } 
    else {
        fDetectionType = BOX_DETECTION;
    }

}

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
    
    G4String volumeName = preStepPoint->GetTouchable()->GetVolume()->GetName();
    
    if (fEventAction->GetRecordAllEvents()) fEventAction->MarkScattering(eventNum);
    
    if (fDetectionType == RPC_DETECTION) {
        if (preStepPoint->GetStepStatus() != fGeomBoundary) return true;
        FillHitData(preStepPoint, eventNum);
    }   
    else if (fDetectionType == BOX_DETECTION) {
        G4bool postStepOnBoundary = (postStepPoint->GetStepStatus() == fGeomBoundary);
        FillScatteringData(preStepPoint, eventNum);
        if (postStepOnBoundary) {
            FillScatteringData(postStepPoint, eventNum);
        }              
    }
    else if (fDetectionType == WORLD_DETECTION) {
        G4bool postStepOnBoundary = (postStepPoint->GetStepStatus() == fGeomBoundary);
        FillScatteringData(preStepPoint, eventNum);
        if (postStepOnBoundary) {
            FillScatteringData(postStepPoint, eventNum);
        }   
    }

 
    return true;
}


void SensitiveDetector::FillHitData(G4StepPoint *stepPoint, G4int eventNum) {
    if (fEventAction->GetRecordHits()){
        HitData data;
        data.eventNum = eventNum;
        data.muonEnergy = stepPoint->GetTotalEnergy() / GeV;
        data.hitTime = stepPoint->GetGlobalTime();
        data.posMuon = stepPoint->GetPosition();
        data.posDetector = stepPoint->GetTouchable()->GetVolume()->GetTranslation();
        data.momMuon = stepPoint->GetMomentum() / GeV;
        fEventAction->CacheRPCHit(eventNum, data);
    }
}


void SensitiveDetector::FillScatteringData(G4StepPoint *stepPoint, G4int eventNum) {
    
    if (fEventAction->GetRecordScatterings()) {
        ScatteringData data;
        data.eventNum = eventNum;
        data.muonEnergy = stepPoint->GetTotalEnergy() / GeV;
        data.posMuon = stepPoint->GetPosition();
        data.momMuon = stepPoint->GetMomentum() / GeV;
        data.scatTime = stepPoint->GetGlobalTime();
        
        const G4Material* material = stepPoint->GetMaterial();
        data.scatMaterial = material->GetName();
        
        fEventAction->CacheScattering(eventNum, data);    
    }    
}








