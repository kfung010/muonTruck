#include "construction.hh"

RPCConstruction::RPCConstruction() {

	fMessenger = new G4GenericMessenger(this, "/cargo/", "Cargo Construction");
	fMessenger->DeclareProperty("cargoMaterial", cargoMaterial, "Material in cargo");  //lithium, aluminium, copper, lead, ammoniumNitrate

	cargoMaterial = "ammoniumNitrate";

	// Dimensions of the world
    xWorldFull = 50*m;
	yWorldFull = 50*m;
	zWorldFull = 50*m;

	// Dimensions of the cargo
	xTruckFull = 2.5*m;
	yTruckFull = 3*m;
	zTruckFull = 10*m;
	wheelRadius = 0.5*m;

	// Size of one pixel, and the numbers of pixel along the x and z directions
	xPixelFull = 24*cm;
	yPixelFull = 2*mm;
	zPixelFull = 24*cm; 
	numX = 61; //405;   
	numZ = 61; //405;  

	// heights of RPC plates from top/bottom of truck 
	//height = {-3*m, 3*m};  
	height = {-2*m, -1.5*m, -1*m, -0.5*m, 0.5*m, 1*m, 1.5*m, 2*m};

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

	lithium = nist->FindOrBuildMaterial("G4_Li");
	aluminium = nist->FindOrBuildMaterial("G4_Al");
    copper = nist->FindOrBuildMaterial("G4_Cu");
    lead = nist->FindOrBuildMaterial("G4_Pb");
}

G4VPhysicalVolume *RPCConstruction::Construct() {

	// Create world
	solidWorld = new G4Box("solidWorld", xWorldFull/2, yWorldFull/2, zWorldFull/2);
	logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);  

	// Create truck
	if (cargoMaterial == "lithium") cargoMat = lithium;
	else if (cargoMaterial == "aluminium") cargoMat = aluminium;
	else if (cargoMaterial == "copper") cargoMat = copper;
	else if (cargoMaterial == "lead") cargoMat = lead;
	else if (cargoMaterial == "ammoniumNitrate") cargoMat = ammoniumNitrate;
	else throw std::runtime_error("Invalid cargo material. Program ended with an error.");

	solidTruck = new G4Box("solidTruck", xTruckFull/2, yTruckFull/2, zTruckFull/2);
	logicTruck = new G4LogicalVolume(solidTruck, cargoMat, "logicTruck");
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


