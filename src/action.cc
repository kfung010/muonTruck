#include "action.hh"


Action::Action() {}
Action::~Action() {}

void Action::BuildForMaster() const {
	CreateNtuple *createNtuple = new CreateNtuple(); 
	SetUserAction(createNtuple);
}

void Action::Build() const {
	muonGenerator *generator = new muonGenerator();
	SetUserAction(generator);
	CreateNtuple *createNtuple = new CreateNtuple(); 
	SetUserAction(createNtuple);

}