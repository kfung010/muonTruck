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
#include "G4GenericMessenger.hh"
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
		G4double getzRPCFull() const { return numZ*zPixelFull; }
		G4double getyHighestRPC() const { return yTruckFull/2 + (*std::max_element(height.begin(), height.end())); }

	private:
		G4double xWorldFull, yWorldFull, zWorldFull;
		G4double xTruckFull, yTruckFull, zTruckFull;
		G4double wheelRadius;
		G4double xPixelFull, yPixelFull, zPixelFull;
		std::vector<G4double> height;
		G4int numX, numZ;

		G4Tubs *solidWheel;
		G4LogicalVolume *logicWheel;
		G4VPhysicalVolume *physWheel;

		G4Box *solidWorld, *solidTruck, *solidPixel;
		G4LogicalVolume *logicWorld, *logicTruck, *logicPixel;
		G4VPhysicalVolume *physWorld, *physTruck, *physPixel;

		void DefineMaterials();
		G4Material *worldMat, *ammoniumNitrate;

		virtual void ConstructSDandField();
};
#endif
