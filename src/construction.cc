#include "construction.hh"

RPCConstruction::RPCConstruction() {

    fMessenger = new G4GenericMessenger(this, "/apparatus/", "Construction of the apparatus");
    fMessenger->DeclareProperty("worldMaterial", worldMaterial, "World material");  //air, vacuum
    fMessenger->DeclareProperty("rpcMaterial", rpcMaterial, "RPC material");  //air, vacuum
    fMessenger->DeclareProperty("cargoMaterial", cargoMaterial, "Material in cargo");  //lithium, aluminium, copper, lead, ammoniumNitrate, carbon, tungsten, water
    fMessenger->DeclareProperty("cargoThickness", cargoThickness, "Cargo thickness in cm");
    fMessenger->DeclareProperty("cargoLength", cargoLength, "Cargo length in cm");
    fMessenger->DeclareProperty("cargoWidth", cargoWidth, "Cargo width in cm");
    
    fMessenger->DeclareProperty("pixelThickness", pixelThickness, "RPC pixel thickness in mm");
    fMessenger->DeclareProperty("pixelLength", pixelLength, "RPC pixel length in cm");
    fMessenger->DeclareProperty("pixelWidth", pixelWidth, "RPC pixel width in cm");
    fMessenger->DeclareProperty("pixelNum1", numX, "Number of RPC pixels in the x direction");
    fMessenger->DeclareProperty("pixelNum2", numY, "Number of RPC pixels in the y direction");
    
    fMessenger->DeclareProperty("stepLimit", stepLimit, "Step limit in mm"); //testing
    
    
    worldMaterial = "air";
    cargoMaterial = "tungsten";
    rpcMaterial = "air";
    cargoThickness = 100;
    cargoLength = 100;
    cargoWidth = 100;
    pixelThickness = 1;
    pixelLength = 995;
    pixelWidth = 995;
    numX = 1;
    numY = 1;
    
    stepLimit = 0.01;

	  // Dimensions of the world
    xWorldFull = 50*m;
    yWorldFull = 50*m;
    zWorldFull = 50*m;

    // heights of RPC plates from top/bottom of truck 
    //height = {-3*m, 3*m};  
    //height = {-2*m, -1.5*m, -1*m, -0.5*m, 0.5*m, 1*m, 1.5*m, 2*m};
    height = {-90*cm, -70*cm, -50*cm, -30*cm, 30*cm, 50*cm, 70*cm, 90*cm};

}

RPCConstruction::~RPCConstruction() {}

void RPCConstruction::DefineMaterials() {

    G4NistManager *nist = G4NistManager::Instance();
    
    
    if (!G4Material::GetMaterial("H2O")) {
        H2O = new G4Material("H2O", 1.000 * g / cm3, 2);
        H2O->AddElement(nist->FindOrBuildElement("H"), 2);
        H2O->AddElement(nist->FindOrBuildElement("O"), 1);
    } else {
        H2O = G4Material::GetMaterial("H2O");
    }

    if (!G4Material::GetMaterial("ammoniumNitrate")) {
        ammoniumNitrate = new G4Material("ammoniumNitrate", 1.72*g/cm3, 3);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("N"), 2);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("H"), 4);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("O"), 3);
    } else {
        ammoniumNitrate = G4Material::GetMaterial("ammoniumNitrate");
    }

    lithium = nist->FindOrBuildMaterial("G4_Li");
    aluminium = nist->FindOrBuildMaterial("G4_Al");
    copper = nist->FindOrBuildMaterial("G4_Cu");
    tungsten = nist->FindOrBuildMaterial("G4_W");
    lead = nist->FindOrBuildMaterial("G4_Pb");
    carbon = nist->FindOrBuildMaterial("G4_C");
    air = nist->FindOrBuildMaterial("G4_AIR");
    vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    
    if (worldMaterial == "air") worldMat = air;
    else if (worldMaterial == "vacuum") worldMat = vacuum;
    else throw std::runtime_error("Invalid world material. Program ended with an error.");
    
    if (rpcMaterial == "air") rpcMat = air;
    else if (rpcMaterial == "vacuum") rpcMat = vacuum;
    else throw std::runtime_error("Invalid RPC material. Program ended with an error.");

    if (cargoMaterial == "lithium") cargoMat = lithium;
    else if (cargoMaterial == "aluminium") cargoMat = aluminium;
    else if (cargoMaterial == "copper") cargoMat = copper;
    else if (cargoMaterial == "tungsten") cargoMat = tungsten;
    else if (cargoMaterial == "lead") cargoMat = lead;
    else if (cargoMaterial == "ammoniumNitrate") cargoMat = ammoniumNitrate;
    else if (cargoMaterial == "carbon") cargoMat = carbon;
    else if (cargoMaterial == "water") cargoMat = H2O;
    else if (cargoMaterial == "air") cargoMat = aluminium;
    else throw std::runtime_error("Invalid cargo material. Program ended with an error.");
}

G4VPhysicalVolume *RPCConstruction::Construct() {

    DefineMaterials();

    // Create world
    solidWorld = new G4Box("solidWorld", xWorldFull/2, yWorldFull/2, zWorldFull/2);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);  

    // Create truck
    xTruckFull = cargoLength*cm;
    yTruckFull = cargoWidth*cm;
    zTruckFull = cargoThickness*cm;     //thickness
    solidTruck = new G4Box("solidTruck", xTruckFull/2, yTruckFull/2, zTruckFull/2);
    logicTruck = new G4LogicalVolume(solidTruck, cargoMat, "logicTruck");
    
    // Set step limits to ensure accurate simulation
    G4UserLimits* userLimits = new G4UserLimits();  
    userLimits->SetMaxAllowedStep(stepLimit*mm);
    logicTruck->SetUserLimits(userLimits);
    
    physTruck = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicTruck, "logicTruck", logicWorld, false, 0, true);

    // Create boxes
//    xBox1Full = cargoLength*cm;
//    zBox1Full = cargoWidth*cm;
//    yBox1Full = cargoThickness*cm;
//    solidBox1 = new G4Box("solidBox1", xBox1Full/2, yBox1Full/2, zBox1Full/2);
//    logicBox1 = new G4LogicalVolume(solidBox1, cargoMat, "logicBox1");
//    
//    xBox2Full = cargoLength*cm;
//    zBox2Full = cargoWidth*cm;
//    yBox2Full = cargoThickness*cm;
//    solidBox2 = new G4Box("solidBox2", xBox2Full/2, yBox2Full/2, zBox2Full/2);
//    logicBox2 = new G4LogicalVolume(solidBox2, cargoMat, "logicBox2");
//    
//    xBox3Full = cargoLength*cm;
//    zBox3Full = cargoWidth*cm;
//    yBox3Full = cargoThickness*cm;
//    solidBox3 = new G4Box("solidBox3", xBox3Full/2, yBox3Full/2, zBox3Full/2);
//    logicBox3 = new G4LogicalVolume(solidBox3, cargoMat, "logicBox3");
    


    // Create RPC "pixels"
    xPixelFull = pixelLength*cm;
    yPixelFull = pixelWidth*cm;
    zPixelFull = pixelThickness*mm;
    solidPixel = new G4Box("solidPixel", xPixelFull/2, yPixelFull/2, zPixelFull/2);
    logicPixel = new G4LogicalVolume(solidPixel, rpcMat, "logicPixel");
    
    for (G4int i = 0; i < numY; i++) {
        for (G4int j = 0; j < numX; j++) {
            G4double xTrans = -(numX-1)/2.*xPixelFull + j * xPixelFull;
            G4double yTrans = -(numY-1)/2.*yPixelFull + i * yPixelFull;
            for (G4int k = 0; k < height.size(); k++) {
                G4double zTrans = (height[k] > 0 ? zTruckFull/2+height[k] : -zTruckFull/2+height[k]);
                //G4double yTrans = (height[k] > 0 ? yBox1Full/2+height[k] : -yBox1Full/2+height[k]);
                G4Translate3D transRPC(G4ThreeVector(xTrans, yTrans, zTrans));
                physPixel = new G4PVPlacement(transRPC, logicPixel, "physPixel", logicWorld, false, numY*numX*k+numX*i+j, false);
            }
        }
    }

    return physWorld;
}

void RPCConstruction::ConstructSDandField() {

    EventAction* eventAction = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());

    SensitiveDetector *rpcdet = new SensitiveDetector("rpcdet", "rpcdet", eventAction);
    logicPixel->SetSensitiveDetector(rpcdet);

    SensitiveDetector *truckdet = new SensitiveDetector("truckdet", "truckdet", eventAction);
    logicTruck->SetSensitiveDetector(truckdet);   
    
    //SensitiveDetector *airdet = new SensitiveDetector("airdet", "airdet", eventAction);
    //logicWorld->SetSensitiveDetector(airdet);   
    
//    SensitiveDetector *box1det = new SensitiveDetector("box1det", "box1det");  
//    logicBox1->SetSensitiveDetector(box1det);   
//    SensitiveDetector *box2det = new SensitiveDetector("box2det", "box2det");  
//    logicBox2->SetSensitiveDetector(box2det);   
//    SensitiveDetector *box3det = new SensitiveDetector("box3det", "box3det");  
//    logicBox3->SetSensitiveDetector(box3det);   
}


