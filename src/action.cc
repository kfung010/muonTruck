#include "action.hh"


Action::Action() {}
Action::~Action() {}

void Action::Build() const {
	muonGenerator *generator = new muonGenerator();
	SetUserAction(generator);
	CreateNtuple *createNtuple = new CreateNtuple(); 
	SetUserAction(createNtuple);

}