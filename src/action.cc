#include "action.hh"


Action::Action() {}
Action::~Action() {}

void Action::BuildForMaster() const {
    CreateNtuple *createNtuple = new CreateNtuple(); 
    SetUserAction(createNtuple);
}

void Action::Build() const {

    EventAction *eventAction = new EventAction();
    
    muonGenerator *generator = new muonGenerator(eventAction);
    SetUserAction(generator);
    SetUserAction(eventAction);
    
    CreateNtuple *createNtuple = new CreateNtuple(); 
    SetUserAction(createNtuple);
    
}