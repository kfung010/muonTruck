#include "eventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"

EventAction::EventAction() {}
EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    fRPCHitCache[eventID].clear();
    fScatteringCache[eventID].clear();
    const CreateNtuple* ntuple = dynamic_cast<const CreateNtuple*>(G4RunManager::GetRunManager()->GetUserRunAction());
    SetRecordAllEvents(ntuple->GetRecordAllEvents());
    SetRecordGenerator(ntuple->GetRecordGenerator());
    SetRecordHits(ntuple->GetRecordHits());
    SetRecordScatterings(ntuple->GetRecordScatterings());
} 
    
    
void EventAction::EndOfEventAction(const G4Event *event) {

    G4int eventID = event->GetEventID();
    if (eventID % 1000 == 0) {
        G4cout << eventID << " events generated ..." << G4endl;
    }
    
    const CreateNtuple* ntuple = dynamic_cast<const CreateNtuple*>(G4RunManager::GetRunManager()->GetUserRunAction());
    
    G4int genId = ntuple->GetGeneratorNtupleId();
    G4int hitsId = ntuple->GetHitsNtupleId();
    G4int scatId = ntuple->GetScatteringsNtupleId();
    
    if (fEventsWithScatterings.count(eventID)) {
        G4AnalysisManager *manager = G4AnalysisManager::Instance();
        if (ntuple->GetRecordGenerator() && genId >= 0) {
            for (const auto& gen : fGeneratorCache[eventID]) {
                manager->FillNtupleIColumn(genId, 0, gen.eventNum);
                manager->FillNtupleDColumn(genId, 1, gen.muonEnergy);
                manager->FillNtupleDColumn(genId, 2, gen.momMuon.x());
                manager->FillNtupleDColumn(genId, 3, gen.momMuon.y());
                manager->FillNtupleDColumn(genId, 4, gen.momMuon.z());
                manager->AddNtupleRow(genId);
            }
        }
        if (ntuple->GetRecordHits() && hitsId >= 0) {
            for (const auto& hit : fRPCHitCache[eventID]) {
                manager->FillNtupleIColumn(hitsId, 0, hit.eventNum);
                manager->FillNtupleDColumn(hitsId, 1, hit.muonEnergy);
                manager->FillNtupleDColumn(hitsId, 2, hit.hitTime);
                manager->FillNtupleDColumn(hitsId, 3, hit.posMuon.x());
                manager->FillNtupleDColumn(hitsId, 4, hit.posMuon.y());
          			manager->FillNtupleDColumn(hitsId, 5, hit.posMuon.z());
          			manager->FillNtupleDColumn(hitsId, 6, hit.posDetector.x());
          			manager->FillNtupleDColumn(hitsId, 7, hit.posDetector.y());
          			manager->FillNtupleDColumn(hitsId, 8, hit.posDetector.z());
          			manager->FillNtupleDColumn(hitsId, 9, hit.momMuon.x());
          			manager->FillNtupleDColumn(hitsId, 10, hit.momMuon.y());
          			manager->FillNtupleDColumn(hitsId, 11, hit.momMuon.z());
                manager->AddNtupleRow(hitsId);
            }
        }
        if (recordScatterings && ntuple->GetRecordScatterings() && scatId >= 0) {
            for (const auto& scat : fScatteringCache[eventID]) {
                manager->FillNtupleIColumn(scatId, 0, scat.eventNum);
                manager->FillNtupleDColumn(scatId, 1, scat.muonEnergy);
                manager->FillNtupleDColumn(scatId, 2, scat.scatTime);
                manager->FillNtupleDColumn(scatId, 3, scat.posMuon.x());
                manager->FillNtupleDColumn(scatId, 4, scat.posMuon.y());
                manager->FillNtupleDColumn(scatId, 5, scat.posMuon.z());
                manager->FillNtupleDColumn(scatId, 6, scat.momMuon.x());
                manager->FillNtupleDColumn(scatId, 7, scat.momMuon.y());
                manager->FillNtupleDColumn(scatId, 8, scat.momMuon.z());
                manager->FillNtupleSColumn(scatId, 9, scat.scatMaterial);
                manager->AddNtupleRow(scatId);
            }
        }
    }
    
    /*G4cout << "Event " << eventID 
       << ": RPC hits = " << fRPCHitCache[eventID].size()
       << ", Scatterings = " << fEventsWithScatterings.count(eventID)
       << G4endl;*/

    fRPCHitCache.erase(eventID);
    fGeneratorCache.erase(eventID);
    fScatteringCache.erase(eventID);
    fEventsWithScatterings.erase(eventID);
}


void EventAction::CacheGenerator(G4int eventNum, const GeneratorData &data) { 
    fGeneratorCache[eventNum].push_back(data);
}
void EventAction::CacheRPCHit(G4int eventNum, const HitData &data) {
    fRPCHitCache[eventNum].push_back(data);
}
void EventAction::CacheScattering(G4int eventNum, const ScatteringData &data) {
    fScatteringCache[eventNum].push_back(data);
}

void EventAction::MarkScattering(G4int eventNum) {
    fEventsWithScatterings.insert(eventNum);
}