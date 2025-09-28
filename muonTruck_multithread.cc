#ifdef _WIN32
  #include <direct.h> // For _mkdir on Windows
  #define mkdir(dir, mode) _mkdir(dir) // ignore mode on Windows
#else
  #include <sys/stat.h>
  #include <sys/types.h>
#endif
#include <string>
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

std::string getOutputDirFromConfig(const char* configFileName) {
    std::string filename(configFileName);
    size_t lastSlash = filename.find_last_of("/\\");
    std::string baseName = (lastSlash == std::string::npos) ? filename : filename.substr(lastSlash + 1);
    size_t dotPos = baseName.find_last_of('.');
    return (dotPos == std::string::npos) ? baseName : baseName.substr(0, dotPos);
}

int main(int argc, char** argv) {


  	#ifdef G4MULTITHREADED
  	    G4MTRunManager* runManager = new G4MTRunManager();
  	#else
  	    G4RunManager* runManager = new G4RunManager();
  	#endif
  	
  
  	runManager->SetUserInitialization(new RPCConstruction());   
  	runManager->SetUserInitialization(new PhysicsList());
  
  	if (argc == 1) {
  		G4cerr << "Please use the configuration file, e.g. runMultiMuon.mac" << G4endl;
          delete runManager;
          return EXIT_FAILURE;
  	}
    
    std::string filepath = argv[1];
    size_t lastdot = filepath.find_last_of(".");
    std::string foldername = (lastdot == std::string::npos) ? filepath : filepath.substr(0, lastdot);
    mkdir(foldername.c_str(), 0755);
    runManager->SetUserInitialization(new Action(foldername));
    
  	G4VisManager* visManager = new G4VisExecutive();
  	visManager->Initialize();
   
  	G4UImanager* UImanager = G4UImanager::GetUIpointer();
  	G4String command = "/control/execute ";
  	G4String fileName = argv[1];
  	UImanager->ApplyCommand(command + fileName);

	std::string baseMacroName = getOutputDirFromConfig(fileName.c_str()) + ".mac";
	std::string newPath = foldername + "/" + baseMacroName;
	if (std::rename(fileName.c_str(), newPath.c_str()) != 0) {
		G4cerr << "Warning: failed to move the configuration file to output folder" << G4endl;
	}

    return 0;
   
}
