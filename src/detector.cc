#include "detector.hh"

RPCDetector::RPCDetector(G4String name) : G4VSensitiveDetector(name) {}

RPCDetector::~RPCDetector(){}

G4bool RPCDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROIhist) {

	// Extract the information of the track segment when it just enters an RPC pixel
	G4Track *track = aStep->GetTrack();

	G4ParticleDefinition* particleDef = track->GetDefinition();
	if (!(particleDef->GetParticleName() == "mu+" || particleDef->GetParticleName() == "mu-")) return true;

	G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

	if (preStepPoint->GetStepStatus() != fGeomBoundary) return true;

	G4int eventNum = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	G4ThreeVector posMuon = preStepPoint->GetPosition();
	G4ThreeVector momMuon = preStepPoint->GetMomentum() / GeV;
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4int copyNo = touchable->GetCopyNumber();
	G4VPhysicalVolume *physVol = touchable->GetVolume();
	G4ThreeVector posDetector = physVol->GetTranslation();
	G4double hitTime = preStepPoint->GetGlobalTime();
	G4double muonMomentum = momMuon.mag();

	//G4cout << "Muon hit time : " << hitTime << G4endl;
	//G4cout << "Muon hit position : (" << posMuon[0] << " mm, " << posMuon[1]/1000 << " m, " << posMuon[2] << " mm)" << G4endl;
	//G4cout << "Muon pixel coordinates (: " << posDetector[0] << " mm, " << posDetector[1]/1000 << " m, " << posDetector[2] << " mm)" << G4endl;
	//G4cout << "Muon momentum : " << muonMomentum << " GeV" << G4endl << G4endl;


	// Fill ntuple branches
	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	manager->FillNtupleIColumn(1, 0, eventNum);
	manager->FillNtupleDColumn(1, 1, muonMomentum);
	manager->FillNtupleDColumn(1, 2, hitTime);
	manager->FillNtupleDColumn(1, 3, posMuon[0]);
	manager->FillNtupleDColumn(1, 4, posMuon[1]);
	manager->FillNtupleDColumn(1, 5, posMuon[2]);
	manager->FillNtupleDColumn(1, 6, posDetector[0]);
	manager->FillNtupleDColumn(1, 7, posDetector[1]);
	manager->FillNtupleDColumn(1, 8, posDetector[2]);
	manager->FillNtupleDColumn(1, 9, momMuon[0]);
	manager->FillNtupleDColumn(1, 10, momMuon[1]);
	manager->FillNtupleDColumn(1, 11, momMuon[2]);
	manager->AddNtupleRow(1);

	return true;
}