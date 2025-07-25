#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"
#include "generator.hh"
#include "createNtuple.hh"
#include "eventAction.hh"

class Action : public G4VUserActionInitialization {
    public:
        Action();
        ~Action();

        virtual void Build() const;
        virtual void BuildForMaster() const;
    private:
        CreateNtuple* fNtuple;  // RunAction
        EventAction* fEventAction;  // EventAction
        muonGenerator* fPrimaryGenerator;  // PrimaryGeneratorAction
};

#endif