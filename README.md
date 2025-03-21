GEANT4 code for muon tomography project
===========================

# Introduction

The repository is to perform GEANT4 simulation for the muon tomography project of the Chinese University of Hong Kong.
 
Several plates of resistive plate chambers (RPC) are placed above and below a cargo filled with high density material. When the cosmic muons pass the cargo, they are scattered or absorbed. By reconstructing the cosmic muon tracks from the RPC hits, we can deduce the material inside the cargo from the scattering angle. 


# Installation of GEANT4

## LINUX

* Pre-installation: 

```sh
sudo apt install cmake cmake-curses-gui gcc g++ libexpat1-dev libxmu-dev libmotif-dev qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
```

* After that, please refer to the instructions in this Youtube tutorial video : https://www.youtube.com/watch?v=4DTumUo3IKw&t=2s 

## Windows

* Pre-requisite: Visual Studio 20XX

* Refer to GEANT4_Windows_Installation.pdf in this repository to install GEANT4

# Running instructions

Clone the repository by

```sh
git clone https://github.com/kfung010/muonTruck.git
```
## LINUX
```sh
cd /path/to/muonTruck
mkdir build
cd build
cmake ..  # Need to configure+build again if you change the source files or add new files
make
./muonTruck    # Or ./muonTruck_multithread runMultiMuon.mac to run in batch mode and multithread mode
```

## Windows
```sh
cd \path\to\muonTruck
mkdir build
cd build
cmake ..  # Need to configure+build again if you change the source files or add new files
msbuild muTruck.sln /p:Configuration=Release /m:4
Release\muonTruck.exe  # Or Release\muonTruck_multithread.exe runMultiMuon.mac to run in batch mode
```

The apparatus and the muon generation can be controlled by the following self-defined macro commands:
- /apparatus/worldMaterial : air, vacuum
- /apparatus/cargoMaterial : lithium, aluminium, copper, lead, ammoniumNitrate, carbon, tungsten, water
- /apparatus/rpcMaterial : air, vacuum
- /apparatus/cargoThickness (in cm)
- /generation/distribution : monoEnergy_vertical, randomEnergy_vertical, randomEnergy_randomAngle, cosmic 

# Documentation
## Folders and files
* `\include` contains the necessary header files, in which the declarations of variables, functions and classes are put
* `\src` contains the source files, in which the definitions and actual implementations of the functions and classes are put
* `muonTruck.cc` or  `muonTruck_multithread_.cc`  are the main files to execute
* `vis.mac` contains the macro commands for the visualization settings in the interaction interface
* `runMultiMuon.mac` contains the macro commands to run in batch mode (i.e. run multiple events and do not show the interaction interface). 

## Modules
* `construction` defines the geometry of the cargo and RPCs, their positions and the materials filling them
* `generator` defines how a muon is generated (e.g. position, momentum, direction)
* `detector` defines which components to detect the muons and extract the track information from those components
* `createNtuple` creates a ROOT output file and store the necessary information during the simulation
* `physics` defines the physics lists, i.e. what physical interactions (e.g. EM interactions, radioactive decays) should be considered in the simulation 
* `action` defines the sequence of actions to be carried out in the simulation, e.g. muon generation, ntuple creation and etc.