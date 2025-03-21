#include "construction.hh"

RPCConstruction::RPCConstruction() {

    fMessenger = new G4GenericMessenger(this, "/apparatus/", "Construction of the apparatus");
    fMessenger->DeclareProperty("worldMaterial", worldMaterial, "World material");  //air, vacuum
    fMessenger->DeclareProperty("rpcMaterial", rpcMaterial, "RPC material");  //air, vacuum
    fMessenger->DeclareProperty("cargoMaterial", cargoMaterial, "Material in cargo");  //lithium, aluminium, copper, lead, ammoniumNitrate, carbon, tungsten, water
    fMessenger->DeclareProperty("cargoThickness", cargoThickness, "Cargo thickness in cm");
    
    worldMaterial = "air";
    cargoMaterial = "water";
    rpcMaterial = "air";
    cargoThickness = 10;

	  // Dimensions of the world
    xWorldFull = 10*m;
    yWorldFull = 10*m;
    zWorldFull = 10*m;

    // Dimensions of the cargo
    xTruckFull = 30*cm;
    zTruckFull = 30*cm;

    // Size of one pixel, and the numbers of pixel along the x and z directions
    xPixelFull = 1*cm;
    yPixelFull = 1*mm;
    zPixelFull = 1*cm; 
    numX = 155;
    numZ = 155;

    // heights of RPC plates from top/bottom of truck 
    //height = {-3*m, 3*m};  
    //height = {-2*m, -1.5*m, -1*m, -0.5*m, 0.5*m, 1*m, 1.5*m, 2*m};
    height = {-90*cm, -70*cm, -50*cm, -30*cm, 30*cm, 50*cm, 70*cm, 90*cm};

}

RPCConstruction::~RPCConstruction() {}

void RPCConstruction::DefineMaterials() {

    G4NistManager *nist = G4NistManager::Instance();
    
    if (worldMaterial == "air") worldMat = nist->FindOrBuildMaterial("G4_AIR");
    else if (worldMaterial == "vacuum") worldMat = nist->FindOrBuildMaterial("G4_Galactic");
    else throw std::runtime_error("Invalid world material. Program ended with an error.");
    
    if (rpcMaterial == "air") rpcMat = nist->FindOrBuildMaterial("G4_AIR");
    else if (rpcMaterial == "vacuum") rpcMat = nist->FindOrBuildMaterial("G4_Galactic");
    else throw std::runtime_error("Invalid RPC material. Program ended with an error.");

    ammoniumNitrate = new G4Material("ammoniumNitrate", 1.72*g/cm3, 3);
    ammoniumNitrate->AddElement(nist->FindOrBuildElement("N"), 2);
    ammoniumNitrate->AddElement(nist->FindOrBuildElement("H"), 4);
    ammoniumNitrate->AddElement(nist->FindOrBuildElement("O"), 3);
 
    H2O = new G4Material("H2O", 1.000 * g / cm3, 2);
    H2O->AddElement(nist->FindOrBuildElement("H"), 2);
    H2O->AddElement(nist->FindOrBuildElement("O"), 1);


    lithium = nist->FindOrBuildMaterial("G4_Li");
    aluminium = nist->FindOrBuildMaterial("G4_Al");
    copper = nist->FindOrBuildMaterial("G4_Cu");
    tungsten = nist->FindOrBuildMaterial("G4_W");
    lead = nist->FindOrBuildMaterial("G4_Pb");
    carbon = nist->FindOrBuildMaterial("G4_C");
    
    if (cargoMaterial == "lithium") cargoMat = lithium;
    else if (cargoMaterial == "aluminium") cargoMat = aluminium;
    else if (cargoMaterial == "copper") cargoMat = copper;
    else if (cargoMaterial == "tungsten") cargoMat = tungsten;
    else if (cargoMaterial == "lead") cargoMat = lead;
    else if (cargoMaterial == "ammoniumNitrate") cargoMat = ammoniumNitrate;
    else if (cargoMaterial == "carbon") cargoMat = carbon;
    else if (cargoMaterial == "water") cargoMat = H2O;
    else throw std::runtime_error("Invalid cargo material. Program ended with an error.");
}

G4VPhysicalVolume *RPCConstruction::Construct() {

    DefineMaterials();

    // Create world
    solidWorld = new G4Box("solidWorld", xWorldFull/2, yWorldFull/2, zWorldFull/2);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);  

    // Create truck
    yTruckFull = cargoThickness*cm;     //thickness
    solidTruck = new G4Box("solidTruck", xTruckFull/2, yTruckFull/2, zTruckFull/2);
    logicTruck = new G4LogicalVolume(solidTruck, cargoMat, "logicTruck");
    
    G4UserLimits* userLimits = new G4UserLimits();
    userLimits->SetMaxAllowedStep(1*m);
    logicTruck->SetUserLimits(userLimits); 
    
    physTruck = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTruck, "physTruck", logicWorld, false, 0, true);

    // Create RPC "pixels"
    solidPixel = new G4Box("solidPixel", xPixelFull/2, yPixelFull/2, zPixelFull/2);
    logicPixel = new G4LogicalVolume(solidPixel, rpcMat, "logicPixel");
    
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
    SensitiveDetector *rpcdet = new SensitiveDetector("rpcdet", "rpcdet");
    logicPixel->SetSensitiveDetector(rpcdet);

    SensitiveDetector *truckdet = new SensitiveDetector("truckdet", "truckdet");
    logicTruck->SetSensitiveDetector(truckdet);     
}


