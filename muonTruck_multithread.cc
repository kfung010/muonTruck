#include <iostream>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "construction.hh"  
#include "physics.hh"  
#include "action.hh"
#include "eventAction.hh"


int main(int argc, char** argv) {


	#ifdef G4MULTITHREADED
	    G4MTRunManager* runManager = new G4MTRunManager();
	#else
	    G4RunManager* runManager = new G4RunManager();
	#endif
	

	runManager->SetUserInitialization(new RPCConstruction());   
	runManager->SetUserInitialization(new PhysicsList());
	runManager->SetUserInitialization(new Action());
 
  //EventAction* eventAction = new EventAction();
	//runManager->SetUserAction(eventAction);

	if (argc == 1) {
		G4cerr << "Please use the configuration file, e.g. runMultiMuon.mac" << G4endl;
        delete runManager;
        return EXIT_FAILURE;
	}

	G4VisManager* visManager = new G4VisExecutive();
	visManager->Initialize();
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4String command = "/control/execute ";
	G4String fileName = argv[1];
	UImanager->ApplyCommand(command + fileName);

	return 0;
}
