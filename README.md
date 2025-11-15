GEANT4 code for muon tomography project
===========================

# Introduction

The repository contains a GEANT4-based simulation for the cosmic muon tomography project at the Chinese University of Hong Kong. The simulation models cosmic muons passing through various high-density materials and being detected by resistive plate chambers (RPCs) placed above and below these materials. As the cosmic muons traverse the materials, they undergo scattering or absorption. By analyzing the muon scattering patterns, we aim to identify and reconstruct the structures of the materials.


# Installation of GEANT4

## LINUX

* Pre-installation: 

```
sudo apt install cmake cmake-curses-gui gcc g++ libexpat1-dev libxmu-dev libmotif-dev qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
```

* After that, please refer to the instructions in this Youtube tutorial video : https://www.youtube.com/watch?v=4DTumUo3IKw&t=2s 

## Windows

* Pre-requisite: Visual Studio 20XX

* Refer to GEANT4_Windows_Installation.pdf in this repository to install GEANT4

# Running instructions

Clone the repository by

```
git clone https://github.com/kfung010/muonTruck.git
```

Create and navigate to the `build` directory:

 ```
cd /path/to/muonTruck
mkdir build
cd build
```

Configure project:

```
cmake ..
```

Build executable:

```
# LINUX
make

# Windows
msbuild muTruck.sln /p:Configuration=Release /m:4
```

Go back to the `muonTruck/` directory and create a `run` directory. Run the simulation in this `run` directory. Two simulation modes are provided: visualization mode and batch-processing mode. In visualization mode, the simulation runs interactively, allowing you to visually observe the cosmic muon trajectories and scatterings. Batch-processing mode runs the simulation automatically without graphical output, enabling large-scale data generation and analysis more efficiently.

For visualization mode, please modify `vis.mac` and run the `muonTruck` macro:

| Mode              | Linux Command                    | Windows Command                          |
|-------------------|----------------------------------|-------------------------------------------|
| Visualization     | `../build/muonTruck`                    | `..\build\Release\muonTruck.exe`                   |

For batch-processing mode, please modify `config.mac`, COPY it to the `run/` directory and rename it to a meaningful name, such as `cosmic_leadBox.mac`. Run the `muonTruck_multithread` macro:

| Mode              | Linux Command                    | Windows Command                          |
|-------------------|----------------------------------|-------------------------------------------|
| Multi-threaded batch processing   | `../build/muonTruck_multithread cosmic_leadBox.mac` | `..\build\Release\muonTruck_multithread.exe cosmic_leadBox.mac` |

After the simulation completes, a directory named after the basename of the config file (`cosmic_leadBox/` in this example), will be created. The configuration file and the output ROOT file will be automatically moved to this folder.


# Documentation
## Main executables
- `muonTruck.cc`: Single-threaded execution with visualization support
- `muonTruck_multithread.cc`: Multi-threaded execution for batch processing

## Core classes
1. Construction of the appratus (`construction.hh/cc`): 

Defines the geometry, positions and the material of the world, RPC plates and high-density mateiral

2. Physics List (`physics.hh/cc`):

Implements the necessary electromagnetic physics processes and includes step limiter for muons

3. Sensitive Detector (`detector.hh/cc`):

Detects muon hits in RPCs and the scattering in volumes

4. Action Initialization (action.hh/cc):

It manages the three levels of actions: run action (`createNtuple`), event action (`eventAction`) and primary generator action (`generator`)

5. Run action (`createNtuple.hh/cc`):

A run is a collection of events. This class defines the structure of the output ROOT ntuples

6. Event action (`eventAction.hh/cc`):

One event refers to one primary particle injection and all its secondary processes. This class controls the filling of the ntuple for each event.

7. Muon Generator (`generator.hh/cc`):

Define the generation of the muons, including the implementation the cosmic muon spectrum

## Configuration files
- `vis.mac`: visualization setup
- `config.mac`: setup for batch processing in multi-thread mode
### Key Configuration Parameters

####  Construction Settings (`/apparatus/`)

| Parameter              | Values                                  | Description                         |
|----------------------|-----------------------------------------|-------------------------------------|
| `worldMaterial`      | `air`, `vacuum`                         | World volume material               |
| `rpcMaterial`        | `air`, `vacuum`, `RPCgas`               | RPC filling                     |
| `pixelThickness`     | `[mm]`                                  | RPC pixel thickness                       |
| `pixelLength`        | `[cm]`                                  | RPC pixel length                          |
| `pixelWidth`        | `[cm]`                                  | RPC pixel width                          |
| `pixelNum1`        |                                  | Number of RPC pixels in the x-direction                          |
| `pixelNum2`        |                                  | Number of RPC pixels in the y-direction                          |
| `addHeight`          | `[cm]`                                  | Z-positions of RPC plates             |
| `addBox`             | `name,l,w,h,x,y,z,material`             | Add box-shaped high-density material (support lithium, aluminium, copper, lead, ammoniumNitrate, carbon, tungsten, water, uranium, iron)               |
| `addHollowBox`             | `name,outerL,outerW,outerH,wallThickness,x,y,z,material`             | Add a hollow box |
| `addEllipsoid`       | `name,xSemi,ySemi,zSemi,x,y,z,material` | Add ellipsoid high-density material                 |
| `addHollowSphere`       | `name,innerR,outerRadius,x,y,z,material` | Add a spherical shell              |
| `addCylinder`        | `name,innerR,outerR,height,x,y,z,material` | Add cylindrical high-density material         |
| `cargoStepLimit`        | `[mm]` | The maximum length that a particle can move in a single step to ensure accurate tracking          |

#### Muon Generation (`/generation/`)

| Parameter            | Values                                  | Description                         |
|----------------------|-----------------------------------------|-------------------------------------|
| `distribution`       | `cosmic`, `monoEnergy_vertical`, `randomEnergy_vertical`, `monoEnergy_tilt`, `randomEnergy_tilt` | Muon energy/angle distribution |
| `muonEnergy`         | `[GeV]`                                 | Energy for mono-energy modes        |
| `muonAngle`          | `[deg]`                                 | Zenith angle for tilted modes       |
| `displacement`       | `[cm]`                                  | Horizontal offset from center in vertical/tilted modes      |

#### Data Output (`/ntuple/`)

| Parameter            | Default                                 | Description                         |
|----------------------|-----------------------------------------|-------------------------------------|
| `recordGenerator`    | `true`                                  | Save the parameters of the generated muons        |
| `recordAllEvents`    | `true`                                  | Save all events (true) or only events hitting the high-density materials (false)                   |
| `recordHits`         | `true`                                  | Save RPC hit information            |
| `recordScatterings`  | `false`                                 | Save detailed scattering data (only for testing as it is slow and makes large files) |

# Reconstruction

Go to the directory containing the config file and the simulation output files, such as `run/cosmic_leadBox/` (you should have `cosmic_leadBox.mac` and `cosmic_leadBox.root` in this directory). Run the POCA reconstruction by

```
python ../../recon/poca.py <basename> <argument1> <argument2> ...
```

For example:
```
python ../../recon/poca.py cosmic_leadBox -r -p -pp -v 5
```


`-r` means to run the POCA point calculation and store the results.

`-p` means to run the processing and store the scattering density profile.

`-pp` means to run the post-processing -- to plot the density distributions.

`-v 5` means setting the voxel size to be 5 cm. The default value will be 10 cm if not specified.


To run the likelihood reconstruction, just run

```
python ../../recon/likelihood.py <basename> <argument1> <argument2> ...
```


