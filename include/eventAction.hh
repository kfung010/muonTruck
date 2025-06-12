#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct GeneratorData {
    G4int eventNum;
    G4double muonMomentum;
    G4ThreeVector momMuon;
};

struct HitData {
    G4int eventNum;
    G4double muonMomentum;
    G4double hitTime;
    G4ThreeVector posMuon;
    G4ThreeVector posDetector;
    G4ThreeVector momMuon;
};

struct TruckData {
    G4int eventNum;
    G4double muonMomentum;
    G4double scatTime;
    G4ThreeVector posMuon;
    G4ThreeVector momMuon;
};

class EventAction : public G4UserEventAction {
	public:
		EventAction();
		virtual ~EventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);
    
    void CacheGenerator(G4int eventNum, const GeneratorData &data);
		void CacheRPCHit(G4int eventNum, const HitData &data);
    void CacheTruckHit(G4int eventNum, const TruckData &data);
		void MarkTruckHit(G4int eventNum);

	private:
    std::unordered_map<G4int, std::vector<GeneratorData>> fGeneratorCache; 
		std::unordered_map<G4int, std::vector<HitData>> fRPCHitCache; 
    std::unordered_map<G4int, std::vector<TruckData>> fTruckHitCache; 
		std::unordered_set<G4int> fEventsWithTruckHits;
};

#endif