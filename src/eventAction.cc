#include "eventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"

EventAction::EventAction() {}
EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    fRPCHitCache[eventID].clear();
    fTruckHitCache[eventID].clear();
} 
    
    
void EventAction::EndOfEventAction(const G4Event *event) {
    G4int eventID = event->GetEventID();
    
    //G4cout << eventID << " fEventsWithTruckHits: " << fEventsWithTruckHits.count(eventID) << G4endl;
    
    if (fEventsWithTruckHits.count(eventID)) {
        G4AnalysisManager *manager = G4AnalysisManager::Instance();
        for (const auto& gen : fGeneratorCache[eventID]) {
            manager->FillNtupleIColumn(0, 0, gen.eventNum);
            manager->FillNtupleDColumn(0, 1, gen.muonMomentum);
            manager->FillNtupleDColumn(0, 2, gen.momMuon.x());
            manager->FillNtupleDColumn(0, 3, gen.momMuon.y());
            manager->FillNtupleDColumn(0, 4, gen.momMuon.z());
            manager->AddNtupleRow(0);
        }
        for (const auto& hit : fRPCHitCache[eventID]) {
            manager->FillNtupleIColumn(1, 0, hit.eventNum);
            manager->FillNtupleDColumn(1, 1, hit.muonMomentum);
            manager->FillNtupleDColumn(1, 2, hit.hitTime);
            manager->FillNtupleDColumn(1, 3, hit.posMuon.x());
            manager->FillNtupleDColumn(1, 4, hit.posMuon.y());
      			manager->FillNtupleDColumn(1, 5, hit.posMuon.z());
      			manager->FillNtupleDColumn(1, 6, hit.posDetector.x());
      			manager->FillNtupleDColumn(1, 7, hit.posDetector.y());
      			manager->FillNtupleDColumn(1, 8, hit.posDetector.z());
      			manager->FillNtupleDColumn(1, 9, hit.momMuon.x());
      			manager->FillNtupleDColumn(1, 10, hit.momMuon.y());
      			manager->FillNtupleDColumn(1, 11, hit.momMuon.z());
            manager->AddNtupleRow(1);
        }
        /*for (const auto& truck : fTruckHitCache[eventID]) {
            manager->FillNtupleIColumn(2, 0, truck.eventNum);
            manager->FillNtupleDColumn(2, 1, truck.muonMomentum);
            manager->FillNtupleDColumn(2, 2, truck.scatTime);
            manager->FillNtupleDColumn(2, 3, truck.posMuon.x());
            manager->FillNtupleDColumn(2, 4, truck.posMuon.y());
            manager->FillNtupleDColumn(2, 5, truck.posMuon.z());
            manager->FillNtupleDColumn(2, 6, truck.momMuon.x());
            manager->FillNtupleDColumn(2, 7, truck.momMuon.y());
            manager->FillNtupleDColumn(2, 8, truck.momMuon.z());
            manager->AddNtupleRow(2);
        }*/
    }
    
    /*G4cout << "Event " << eventID 
       << ": RPC hits = " << fRPCHitCache[eventID].size()
       << ", Truck hit = " << fEventsWithTruckHits.count(eventID)
       << G4endl;*/

    fRPCHitCache.erase(eventID);
    fEventsWithTruckHits.erase(eventID);
}

void EventAction::CacheGenerator(G4int eventNum, const GeneratorData &data) {
    fGeneratorCache[eventNum].push_back(data);
}
void EventAction::CacheRPCHit(G4int eventNum, const HitData &data) {
    fRPCHitCache[eventNum].push_back(data);
}
void EventAction::CacheTruckHit(G4int eventNum, const TruckData &data) {
    fTruckHitCache[eventNum].push_back(data);
}

void EventAction::MarkTruckHit(G4int eventNum) {
    fEventsWithTruckHits.insert(eventNum);
}