#ifndef CREATENTUPLE_HH
#define CREATENTUPLE_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4GenericMessenger.hh"

class CreateNtuple : public G4UserRunAction {

    public:
        CreateNtuple();
        ~CreateNtuple();

        virtual void BeginOfRunAction(const G4Run *);
        virtual void EndOfRunAction(const G4Run *);
        
        G4int GetGeneratorNtupleId() const { return generatorNtupleId_; }
        G4int GetHitsNtupleId() const { return hitsNtupleId_; }
        G4int GetScatteringsNtupleId() const { return scatteringsNtupleId_; }
        
        void SetRecordAllEvents(G4bool record) { recordAllEvents = record; }
        void SetRecordGenerator(G4bool record) { recordGenerator = record; }
        void SetRecordHits(G4bool record) { recordHits = record; }
        void SetRecordScatterings(G4bool record) { recordScatterings = record; }
        
        G4bool GetRecordAllEvents() const { return recordAllEvents; }
        G4bool GetRecordGenerator() const { return recordGenerator; }
        G4bool GetRecordHits() const { return recordHits; }
        G4bool GetRecordScatterings() const { return recordScatterings; }
        
        void SetOutputDirectory(const G4String& dir) { outputDir = dir; }


    private:
        G4bool recordAllEvents = true;
        G4bool recordGenerator = true;
        G4bool recordHits = true;
        G4bool recordScatterings = false;
        
        G4int generatorNtupleId_ = -1;
        G4int hitsNtupleId_ = -1;
        G4int scatteringsNtupleId_ = -1;
        
        G4String outputDir;
};


#endif