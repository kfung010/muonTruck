#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "createNtuple.hh"
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct GeneratorData {
    G4int eventNum;
    G4double muonEnergy;
    G4ThreeVector momMuon;
};

struct HitData {
    G4int eventNum;
    G4double muonEnergy;
    G4double hitTime;
    G4ThreeVector posMuon;
    G4ThreeVector posDetector;
    G4ThreeVector momMuon;
};

struct ScatteringData {
    G4int eventNum;
    G4double muonEnergy;
    G4double scatTime;
    G4ThreeVector posMuon;
    G4ThreeVector momMuon;
    G4String scatMaterial;
};

class EventAction : public G4UserEventAction {
	public:
		EventAction();
		virtual ~EventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);
    
    void CacheGenerator(G4int eventNum, const GeneratorData &data);
		void CacheRPCHit(G4int eventNum, const HitData &data);
    void CacheScattering(G4int eventNum, const ScatteringData &data);
		void MarkScattering(G4int eventNum);
    
    void SetRecordAllEvents(G4bool record) { recordAllEvents = record; }
    void SetRecordGenerator(G4bool record) { recordGenerator = record; }
    void SetRecordHits(G4bool record) { recordHits = record; }
    void SetRecordScatterings(G4bool record) { recordScatterings = record; }
    
    G4bool GetRecordAllEvents() const { return recordAllEvents; }
    G4bool GetRecordGenerator() const { return recordGenerator; }
    G4bool GetRecordHits() const { return recordHits; }
    G4bool GetRecordScatterings() const { return recordScatterings; }

	private:
    std::unordered_map<G4int, std::vector<GeneratorData>> fGeneratorCache; 
		std::unordered_map<G4int, std::vector<HitData>> fRPCHitCache; 
    std::unordered_map<G4int, std::vector<ScatteringData>> fScatteringCache; 
		std::unordered_set<G4int> fEventsWithScatterings;
    G4bool recordAllEvents = true;
    G4bool recordGenerator = true;
    G4bool recordHits = true;
    G4bool recordScatterings = false;
};

#endif