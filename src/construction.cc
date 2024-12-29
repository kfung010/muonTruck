#include "construction.hh"

RPCConstruction::RPCConstruction() {

    xWorldFull = 50*m;
	yWorldFull = 50*m;
	zWorldFull = 50*m;

	xTruckFull = 2.5*m;
	yTruckFull = 3*m;
	zTruckFull = 10*m;
	wheelRadius = 0.5*m;

	xPixelFull = 6*cm;
	yPixelFull = 1*cm;
	zPixelFull = 12*cm; 
	numX = 155;   
	numZ = 155;  

	// heights of RPC plates from top/bottom of truck 
	//height = {-3*m, 3*m};  
	height = {-4*m, -2*m, 2*m, 4*m};

	DefineMaterials();
}
RPCConstruction::~RPCConstruction() {}

void RPCConstruction::DefineMaterials() {

	G4NistManager *nist = G4NistManager::Instance();

	worldMat = nist->FindOrBuildMaterial("G4_AIR");

	ammoniumNitrate = new G4Material("ammoniumNitrate", 1.72*g/cm3, 3);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("N"), 2);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("H"), 4);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("O"), 3);

}

G4VPhysicalVolume *RPCConstruction::Construct() {

	// Create world
	solidWorld = new G4Box("solidWorld", xWorldFull/2, yWorldFull/2, zWorldFull/2);
	logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);  

	// Create truck
	solidTruck = new G4Box("solidTruck", xTruckFull/2, yTruckFull/2, zTruckFull/2);
	logicTruck = new G4LogicalVolume(solidTruck, ammoniumNitrate, "logicTruck");
	physTruck = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTruck, "physTruck", logicWorld, false, 0, true); 

	// Create wheel
	/*G4Rotate3D rotY(90*deg, G4ThreeVector(0,1,0));
	G4Translate3D trans1(G4ThreeVector(0., -yTruckFull/2-wheelRadius, zTruckFull/2-1*m));
	G4Translate3D trans2(G4ThreeVector(0., -yTruckFull/2-wheelRadius, -zTruckFull/2+1*m));
	solidWheel = new G4Tubs("solidWheel", 0., wheelRadius, xTruckFull/2, 0*deg, 360*deg);
	logicWheel = new G4LogicalVolume(solidWheel, worldMat, "logicWheel");
	physWheel = new G4PVPlacement((trans1)*(rotY), logicWheel, "physWheel", logicWorld, false, 0, true); 
	physWheel = new G4PVPlacement((trans2)*(rotY), logicWheel, "physWheel", logicWorld, false, 0, true); */


	// Create RPC "pixels"
	solidPixel = new G4Box("solidPixel", xPixelFull/2, yPixelFull/2, zPixelFull/2);
	logicPixel = new G4LogicalVolume(solidPixel, worldMat, "logicPixel");
		
	for (G4int i = 0; i < numZ; i++) {
		for (G4int j = 0; j < numX; j++) {
			G4double xTrans = -(numX-1)/2.*xPixelFull + j * xPixelFull;
			G4double zTrans = -(numZ-1)/2.*zPixelFull + i * zPixelFull;
			for (G4int k = 0; k < height.size(); k++) {
				G4double yTrans = (height[k] > 0 ? yTruckFull/2+height[k] : -yTruckFull/2+height[k]);
				G4Translate3D transRPC(G4ThreeVector(xTrans, yTrans, zTrans));
				physPixel = new G4PVPlacement(transRPC, logicPixel, "physPixel", logicWorld, false, numZ*numX*k+numX*i+j, false);
			}
		}
	}
	
	return physWorld;
}

void RPCConstruction::ConstructSDandField() {

	RPCDetector *rpcdet = new RPCDetector("rpcdet");
	logicPixel->SetSensitiveDetector(rpcdet);

}


