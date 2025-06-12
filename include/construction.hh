#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "detector.hh"
#include "eventAction.hh"
#include "G4GenericMessenger.hh"
#include "G4SDManager.hh"
#include <G4UserLimits.hh>
#include <vector>


class RPCConstruction : public G4VUserDetectorConstruction {
    public:
        RPCConstruction();
        ~RPCConstruction();
        virtual G4VPhysicalVolume *Construct();   
        G4double getxWorldFull() const { return xWorldFull; }
        G4double getyWorldFull() const { return yWorldFull; }
        G4double getzWorldFull() const { return zWorldFull; }
        G4double getxRPCFull() const { return numX*xPixelFull; }
        G4double getyRPCFull() const { return numY*yPixelFull; }
        
        G4double getxTruckFull() const { return xTruckFull; }
        G4double getyTruckFull() const { return yTruckFull; }
		G4double getzHighestTruck() const { return zTruckFull/2; }
        
        G4double getzHighestRPC() const { return zTruckFull/2 + (*std::max_element(height.begin(), height.end())); }

    private:
        G4double xWorldFull, yWorldFull, zWorldFull;
        G4double xTruckFull, yTruckFull, zTruckFull;
        
        G4double xBox1Full, yBox1Full, zBox1Full;
        G4double xBox2Full, yBox2Full, zBox2Full;
        G4double xBox3Full, yBox3Full, zBox3Full;
        
        G4double wheelRadius;
        G4double xPixelFull, yPixelFull, zPixelFull;
        std::vector<G4double> height;
        G4int numX, numY;

        G4Tubs *solidWheel;
        G4LogicalVolume *logicWheel;
        G4VPhysicalVolume *physWheel;

        G4Box *solidWorld, *solidTruck, *solidPixel;
        G4LogicalVolume *logicWorld, *logicTruck, *logicPixel;
        G4VPhysicalVolume *physWorld, *physTruck, *physPixel;
        
        G4Box *solidBox1, *solidBox2, *solidBox3;
        G4LogicalVolume *logicBox1, *logicBox2, *logicBox3;
        G4VPhysicalVolume *physBox1, *physBox2, *physBox3;
        
        G4double cargoThickness, cargoLength, cargoWidth, pixelThickness, pixelLength, pixelWidth, stepLimit;

        void DefineMaterials();
        G4Material *worldMat;
        G4Material *cargoMat;
        G4Material *rpcMat;
        G4String cargoMaterial, worldMaterial, rpcMaterial;
        
        G4Material *ammoniumNitrate, *H2O, *air, *vacuum;
        G4Material *lithium, *aluminium, *copper, *lead, *tungsten, *carbon;

        virtual void ConstructSDandField();
        
        G4GenericMessenger *fMessenger;
};
#endif
