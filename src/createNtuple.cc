#include "createNtuple.hh"

CreateNtuple::CreateNtuple() {
    
    G4GenericMessenger* messenger = new G4GenericMessenger(this, "/ntuple/", "NTuple configuration");
    
    messenger->DeclareProperty("recordAllEvents", recordAllEvents, "Record all events (true) or only events with scattering in cargo (false)");
    messenger->DeclareProperty("recordGenerator", recordGenerator, "Record generator data");
    messenger->DeclareProperty("recordHits", recordHits, "Record RPC hit data");
    messenger->DeclareProperty("recordScatterings", recordScatterings, "Record scattering data");

}

CreateNtuple::~CreateNtuple() {}

void CreateNtuple::BeginOfRunAction(const G4Run *run) {
    
    G4AnalysisManager *manager = G4AnalysisManager::Instance();
    
    G4int generatorNtupleId = -1;
    G4int hitsNtupleId = -1;
    G4int scatteringsNtupleId = -1;
    
    if (recordGenerator) {
        generatorNtupleId = manager->CreateNtuple("generator", "generator");   // Information of generated muons
        manager->CreateNtupleIColumn(generatorNtupleId, "eventNumber");
        manager->CreateNtupleDColumn(generatorNtupleId, "muonEnergy");
        manager->CreateNtupleDColumn(generatorNtupleId, "muonMomentumX");
        manager->CreateNtupleDColumn(generatorNtupleId, "muonMomentumY");
        manager->CreateNtupleDColumn(generatorNtupleId, "muonMomentumZ");
        manager->FinishNtuple(generatorNtupleId);
    }
    
    if (recordHits) {
        hitsNtupleId = manager->CreateNtuple("hits", "hits");   // Information of RPC hits
        manager->CreateNtupleIColumn(hitsNtupleId, "eventNumber");
        manager->CreateNtupleDColumn(hitsNtupleId, "muonEnergy");  // GeV
        manager->CreateNtupleDColumn(hitsNtupleId, "hitTime");  // ns
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPositionX_truth");  // mm
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPositionY_truth");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPositionZ_truth");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPixelX");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPixelY");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitPixelZ");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitMomentumX_truth");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitMomentumY_truth");
        manager->CreateNtupleDColumn(hitsNtupleId, "hitMomentumZ_truth");
        manager->FinishNtuple(hitsNtupleId);
    }
    
    if (recordScatterings) {
        scatteringsNtupleId = manager->CreateNtuple("scatterings", "scatterings");   // Information of scattering points
        manager->CreateNtupleIColumn(scatteringsNtupleId, "eventNumber");
        manager->CreateNtupleDColumn(scatteringsNtupleId, "muonEnergy");  // GeV
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatTime");  // ns
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatPositionX_truth");  // mm
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatPositionY_truth");
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatPositionZ_truth");
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatMomentumX_truth");
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatMomentumY_truth");
        manager->CreateNtupleDColumn(scatteringsNtupleId, "scatMomentumZ_truth");
        manager->CreateNtupleSColumn(scatteringsNtupleId, "scatMaterial");
        manager->FinishNtuple(scatteringsNtupleId);
    }
    
    generatorNtupleId_ = generatorNtupleId;
    hitsNtupleId_ = hitsNtupleId;
    scatteringsNtupleId_ = scatteringsNtupleId;
        
    G4int runID = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runID;

    G4String fileName = "output_"+strRunID.str()+".root";
    if (!outputDir.empty()) {
        if (outputDir.back() != '/') outputDir += '/';
        fileName = outputDir + fileName;
    }
    manager->OpenFile(fileName);
}

void CreateNtuple::EndOfRunAction(const G4Run *) {
    G4AnalysisManager *manager = G4AnalysisManager::Instance();
    manager->Write();
    manager->CloseFile();
}