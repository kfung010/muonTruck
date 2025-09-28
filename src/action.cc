#include "action.hh"


Action::Action(const G4String& outputDir) : fOutputDir(outputDir) {
    fNtuple = new CreateNtuple();
    fNtuple->SetOutputDirectory(fOutputDir);
}

Action::~Action() {
    delete fNtuple;
}

void Action::Build() const {
    CreateNtuple* workerNtuple = new CreateNtuple();
    workerNtuple->SetOutputDirectory(fOutputDir);
    
    EventAction* workerEventAction = new EventAction();
    muonGenerator* workerGenerator = new muonGenerator(workerEventAction);

    //workerEventAction->SetRecordAllEvents(workerNtuple->GetRecordAllEvents());
    //workerEventAction->SetRecordTruckHits(workerNtuple->GetRecordTruckHits());

    SetUserAction(workerNtuple);            // RunAction
    SetUserAction(workerEventAction);       // EventAction
    SetUserAction(workerGenerator);         // PrimaryGeneratorAction
}

void Action::BuildForMaster() const {
    fNtuple->SetOutputDirectory(fOutputDir);
    SetUserAction(fNtuple);
}