GEANT4 code for muon tomography project
===========================

# Introduction

The repository contains a GEANT4-based simulation for the cosmic muon tomography project of the Chinese University of Hong Kong. The simulation models cosmic muons passing through varius high-density materials and being detected by the resistive plate chambers (RPC) placed above and below the materials. When the cosmic muons pass the materials, they are scattered or absorbed. By analyzing muon scattering patterns, we hope to identify and  and reconstruct the 3D structures of the materials.


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

Run the simulation in the `build` directory. See below for the options in the `config.mac` file.

| Mode              | Linux Command                    | Windows Command                          |
|-------------------|----------------------------------|-------------------------------------------|
| Visualization     | `./muonTruck`                    | `Release\muonTruck.exe`                   |
| Multi-threaded batch processing   | `./muonTruck_multithread ../config.mac` | `Release\muonTruck_multithread.exe ../config.mac` |



# Documentation
## Main executables
- `muonTruck.cc`: Single-threaded execution with visualization support
- `muonTruck_multithread.cc`: Multi-threaded execution for batch processing

## Core classes
1. Construction of the appratus (`construction.hh/cc`)
Defines the geometry, positions and the material of the world, RPC plates and high-density mateiral
2. Physics List (`physics.hh/cc`)
Implements the necessary electromagnetic physics processes and includes step limiter for muons
3. Sensitive Detector (`detector.hh/cc`)
Detects muon hits in RPCs and the scattering in volumes
4. Action Initialization (action.hh/cc)
It manages the three levels of actions: run action (`createNtuple`), event action (`eventAction`) and primary generator action (`generator`)
5. Run action (`createNtuple.hh/cc`)
A run is a collection of events. This class defines the structure of the output ROOT ntuples
6. Event action (`eventAction.hh/cc`)
One event refers to one primary particle injection and all its secondary processes. This class controls the filling of the ntuple for each event.
7. Muon Generator (`generator.hh/cc`)
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
| `addBox`             | `name,l,w,h,x,y,z,material`             | Add box-shaped high-density material (support lithium, aluminium, copper, lead, ammoniumNitrate, carbon, tungsten, water)               |
| `addEllipsoid`       | `name,xSemi,ySemi,zSemi,x,y,z,material` | Add ellipsoid high-density material                 |
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