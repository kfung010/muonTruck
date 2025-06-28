#include "action.hh"


Action::Action() {
    fNtuple = new CreateNtuple();
}

Action::~Action() {
    delete fNtuple;
}

void Action::Build() const {
    CreateNtuple* workerNtuple = new CreateNtuple();
    EventAction* workerEventAction = new EventAction();
    muonGenerator* workerGenerator = new muonGenerator(workerEventAction);

    //workerEventAction->SetRecordAllEvents(workerNtuple->GetRecordAllEvents());
    //workerEventAction->SetRecordTruckHits(workerNtuple->GetRecordTruckHits());

    SetUserAction(workerNtuple);            // RunAction
    SetUserAction(workerEventAction);       // EventAction
    SetUserAction(workerGenerator);         // PrimaryGeneratorAction
}

void Action::BuildForMaster() const {
    SetUserAction(fNtuple);
}