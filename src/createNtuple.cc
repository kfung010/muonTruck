#include "createNtuple.hh"

CreateNtuple::CreateNtuple() {

    G4AnalysisManager *manager = G4AnalysisManager::Instance();

    manager->CreateNtuple("generator", "generator");   // Information of generated muons
    manager->CreateNtupleDColumn("muonMomentum");
    manager->CreateNtupleDColumn("muonMomentumX");
    manager->CreateNtupleDColumn("muonMomentumY");
    manager->CreateNtupleDColumn("muonMomentumZ");
    manager->CreateNtupleDColumn("muonMomentum_Zenith");
    manager->FinishNtuple(0);

    manager->CreateNtuple("hits", "hits");   // Information of RPC hits
    manager->CreateNtupleIColumn("eventNumber");
    manager->CreateNtupleDColumn("muonMomentum");  // GeV
    manager->CreateNtupleDColumn("hitTime");  // ns
    manager->CreateNtupleDColumn("hitPositionX_truth");  // mm
    manager->CreateNtupleDColumn("hitPositionY_truth");
    manager->CreateNtupleDColumn("hitPositionZ_truth");
    manager->CreateNtupleDColumn("hitPixelX");
    manager->CreateNtupleDColumn("hitPixelY");
    manager->CreateNtupleDColumn("hitPixelZ");
    manager->CreateNtupleDColumn("hitMomentumX_truth");
    manager->CreateNtupleDColumn("hitMomentumY_truth");
    manager->CreateNtupleDColumn("hitMomentumZ_truth");
    manager->FinishNtuple(1);
    
    manager->CreateNtuple("truck", "truck");   // Information of scattering points in the truck
    manager->CreateNtupleIColumn("eventNumber");
    manager->CreateNtupleDColumn("muonMomentum");  // GeV
    manager->CreateNtupleDColumn("scatTime");  // ns
    manager->CreateNtupleDColumn("scatPositionX_truth");  // mm
    manager->CreateNtupleDColumn("scatPositionY_truth");
    manager->CreateNtupleDColumn("scatPositionZ_truth");
    manager->CreateNtupleDColumn("scatMomentumX_truth");
    manager->CreateNtupleDColumn("scatMomentumY_truth");
    manager->CreateNtupleDColumn("scatMomentumZ_truth");
    manager->FinishNtuple(2);

}

CreateNtuple::~CreateNtuple() {}

void CreateNtuple::BeginOfRunAction(const G4Run *run) {

    G4AnalysisManager *manager = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runID;

    G4String fileName = "output_"+strRunID.str()+".root";
    manager->OpenFile(fileName);
}

void CreateNtuple::EndOfRunAction(const G4Run *) {
    G4AnalysisManager *manager = G4AnalysisManager::Instance();
    manager->Write();
    manager->CloseFile();
}