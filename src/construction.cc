#include "construction.hh"

RPCConstruction::RPCConstruction() {
    
    DefineMaterials();

    fMessenger = new G4GenericMessenger(this, "/apparatus/", "Construction of the apparatus");
    fMessenger->DeclareProperty("worldMaterial", worldMaterial, "World material");
    fMessenger->DeclareProperty("rpcMaterial", rpcMaterial, "RPC material");
    fMessenger->DeclareProperty("pixelThickness", pixelThickness, "RPC pixel thickness in mm");
    fMessenger->DeclareProperty("pixelLength", pixelLength, "RPC pixel length in cm");
    fMessenger->DeclareProperty("pixelWidth", pixelWidth, "RPC pixel width in cm");
    fMessenger->DeclareProperty("pixelNum1", numX, "Number of RPC pixels in the x direction");
    fMessenger->DeclareProperty("pixelNum2", numY, "Number of RPC pixels in the y direction");
    fMessenger->DeclareProperty("cargoStepLimit", cargoStepLimit, "Cargo step limit in mm");
    fMessenger->DeclareProperty("surroundingStepLimit", surroundingStepLimit, "Surroundings step limit in mm");
    
    fMessenger->DeclareMethod("clearHeights", 
                             &RPCConstruction::ClearRPCHeights, 
                             "Clear all RPC heights");
    
    fMessenger->DeclareMethod("addHeight", 
                             &RPCConstruction::AddRPCHeight, 
                             "Add RPC height (in cm)");
    
    fMessenger->DeclareMethod("clearCargo", 
                             &RPCConstruction::ClearCargo, 
                             "Clear all cargo");
    
    fMessenger->DeclareMethod("addBox", 
                             &RPCConstruction::AddCargoBox, 
                             "Add box cargo: name l(cm) w(cm) h(cm) x(cm) y(cm) z(cm) material");
    
    fMessenger->DeclareMethod("addEllipsoid", 
                             &RPCConstruction::AddCargoEllipsoid, 
                             "Add ellipsoid cargo: name xSemiAxis(cm) ySemiAxis(cm) zSemiAxis(cm) x(cm) y(cm) z(cm) material");
    
    fMessenger->DeclareMethod("addCylinder", 
                             &RPCConstruction::AddCargoCylinder, 
                             "Add cylinder cargo: name innerRadius(cm) outerRadius(cm) height(cm) x(cm) y(cm) z(cm) material");
    
}

RPCConstruction::~RPCConstruction() {}

void RPCConstruction::DefineMaterials() {

    G4NistManager *nist = G4NistManager::Instance();
    
    
    if (!G4Material::GetMaterial("H2O", false)) {
        H2O = new G4Material("H2O", 1.000 * g / cm3, 2);
        H2O->AddElement(nist->FindOrBuildElement("H"), 2);
        H2O->AddElement(nist->FindOrBuildElement("O"), 1);
    } 
    else {
        H2O = G4Material::GetMaterial("H2O");
    }

    if (!G4Material::GetMaterial("ammoniumNitrate", false)) {
        ammoniumNitrate = new G4Material("ammoniumNitrate", 1.72*g/cm3, 3);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("N"), 2);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("H"), 4);
        ammoniumNitrate->AddElement(nist->FindOrBuildElement("O"), 3);
    } 
    else {
        ammoniumNitrate = G4Material::GetMaterial("ammoniumNitrate");
    }
    
    if (!G4Material::GetMaterial("C2H2F4", false)) {
        C2H2F4 = new G4Material("C2H2F4", 4.25*mg/cm3, 3, kStateGas, 293.15*kelvin, 1*atmosphere);
        C2H2F4->AddElement(nist->FindOrBuildElement("C"), 2);
        C2H2F4->AddElement(nist->FindOrBuildElement("H"), 2);
        C2H2F4->AddElement(nist->FindOrBuildElement("F"), 4);
    } 
    else {
        C2H2F4 = G4Material::GetMaterial("C2H2F4");
    }
    
    if (!G4Material::GetMaterial("iC4H10", false)) {
        iC4H10 = new G4Material("iC4H10", 2.51*mg/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        iC4H10->AddElement(nist->FindOrBuildElement("C"), 4);
        iC4H10->AddElement(nist->FindOrBuildElement("H"), 10);
    } 
    else {
        iC4H10 = G4Material::GetMaterial("iC4H10");
    }
    
    if (!G4Material::GetMaterial("RPCgas", false)) {
        G4double density = 0.95*4.25*mg/cm3 + 0.05*2.51*mg/cm3;
        RPCgas = new G4Material("RPCgas", density, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        RPCgas->AddMaterial(C2H2F4, 0.95);
        RPCgas->AddMaterial(iC4H10, 0.05);
    } 
    else {
        RPCgas = G4Material::GetMaterial("RPCgas");
    }

    lithium = nist->FindOrBuildMaterial("G4_Li");
    aluminium = nist->FindOrBuildMaterial("G4_Al");
    copper = nist->FindOrBuildMaterial("G4_Cu");
    iron = nist->FindOrBuildMaterial("G4_Fe");
    uranium = nist->FindOrBuildMaterial("G4_U");
    tungsten = nist->FindOrBuildMaterial("G4_W");
    lead = nist->FindOrBuildMaterial("G4_Pb");
    carbon = nist->FindOrBuildMaterial("G4_C");
    air = nist->FindOrBuildMaterial("G4_AIR");
    vacuum = nist->FindOrBuildMaterial("G4_Galactic");

}

void RPCConstruction::ClearRPCHeights() {
    rpcHeights.clear();
}

void RPCConstruction::AddRPCHeight(G4double height) {
    rpcHeights.push_back(height * cm);
}

void RPCConstruction::AddCargoBox(const G4String& name, G4double length, G4double width, 
                                 G4double height, G4double posX, G4double posY, G4double posZ,
                                 const G4String& materialName) {
    CargoBox box;
    box.name = name;
    box.width = width;
    box.height = height;
    box.length = length;
    box.position = G4ThreeVector(posX, posY, posZ);
    box.material = GetMaterialByName(materialName);
    cargoBoxes[name] = box;
}

void RPCConstruction::AddCargoEllipsoid(const G4String& name, 
                                        G4double xSemiAxis, G4double ySemiAxis, G4double zSemiAxis,
                                        G4double posX, G4double posY, G4double posZ,
                                        const G4String& materialName) {
    CargoEllipsoid ellipsoid;
    ellipsoid.name = name;
    ellipsoid.xSemiAxis = xSemiAxis;
    ellipsoid.ySemiAxis = ySemiAxis;
    ellipsoid.zSemiAxis = zSemiAxis;
    ellipsoid.position = G4ThreeVector(posX, posY, posZ);
    ellipsoid.material = GetMaterialByName(materialName);
    cargoEllipsoids[name] = ellipsoid;
}

void RPCConstruction::AddCargoCylinder(const G4String& name, 
                                       G4double innerRadius, G4double outerRadius, G4double height,
                                       G4double posX, G4double posY, G4double posZ,
                                       const G4String& materialName) {
    CargoCylinder cylinder;
    cylinder.name = name;
    cylinder.innerRadius = innerRadius;
    cylinder.outerRadius = outerRadius;
    cylinder.height = height;
    cylinder.position = G4ThreeVector(posX, posY, posZ);
    cylinder.material = GetMaterialByName(materialName);
    cargoCylinders[name] = cylinder;
}


void RPCConstruction::ClearCargo() {
    cargoBoxes.clear();
    cargoEllipsoids.clear();
    cargoCylinders.clear();

}

void RPCConstruction::AddCargoBox(const G4String& params) {
    std::vector<G4String> tokens;
    std::istringstream iss(params);
    G4String token;
    
    while (std::getline(iss, token, ',')) {
        size_t start = token.find_first_not_of(" ");
        size_t end = token.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            tokens.push_back(token.substr(start, end - start + 1));
        } else if (!token.empty()) {
            tokens.push_back(token);
        }
    }

    if (tokens.size() < 8) {
        G4cerr << "Error: Insufficient parameters for addBox. Expected 8, got " 
               << tokens.size() << G4endl;
        return;
    }
    
    try {
        G4double l = std::stod(tokens[1]);
        G4double w = std::stod(tokens[2]);
        G4double h = std::stod(tokens[3]);
        G4double x = std::stod(tokens[4]);
        G4double y = std::stod(tokens[5]);
        G4double z = std::stod(tokens[6]);
        
        G4cout << "AddCargoBox: " << tokens[0] << " " << l << " " << w << " " << h 
               << " " << x << " " << y << " " << z << " " << tokens[7] << G4endl;
        
        AddCargoBox(tokens[0], l*cm, w*cm, h*cm, x*cm, y*cm, z*cm, tokens[7]);
    }
    catch (const std::exception& e) {
        G4cerr << "Error parsing parameters: " << e.what() << G4endl;
    }
}

void RPCConstruction::AddCargoEllipsoid(const G4String& params) {
    std::vector<G4String> tokens;
    std::istringstream iss(params);
    G4String token;
    
    while (std::getline(iss, token, ',')) {
        size_t start = token.find_first_not_of(" ");
        size_t end = token.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            tokens.push_back(token.substr(start, end - start + 1));
        } else if (!token.empty()) {
            tokens.push_back(token);
        }
    }

    if (tokens.size() < 8) {
        G4cerr << "Error: Insufficient parameters for addEllipsoid. Expected 8, got " 
               << tokens.size() << G4endl;
        return;
    }
    
    try {
        G4double xSemi = std::stod(tokens[1]);
        G4double ySemi = std::stod(tokens[2]);
        G4double zSemi = std::stod(tokens[3]);
        G4double x = std::stod(tokens[4]);
        G4double y = std::stod(tokens[5]);
        G4double z = std::stod(tokens[6]);
        
        G4cout << "AddCargoEllipsoid: " << tokens[0] << " " << xSemi << " " << ySemi << " " << zSemi 
               << " " << x << " " << y << " " << z << " " << tokens[7] << G4endl;
        
        AddCargoEllipsoid(tokens[0], 
                         xSemi*cm, ySemi*cm, zSemi*cm,
                         x*cm, y*cm, z*cm,
                         tokens[7]);
    }
    catch (const std::exception& e) {
        G4cerr << "Error parsing parameters: " << e.what() << G4endl;
    }
}

void RPCConstruction::AddCargoCylinder(const G4String& params) {

    std::vector<G4String> tokens;
    std::istringstream iss(params);
    G4String token;
    
    while (std::getline(iss, token, ',')) {
        size_t start = token.find_first_not_of(" ");
        size_t end = token.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            tokens.push_back(token.substr(start, end - start + 1));
        } else if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    
    if (tokens.size() < 8) {
        G4cerr << "Error: Insufficient parameters for addCylinder. Expected 8, got " 
               << tokens.size() << G4endl;
        return;
    }
    
    try {
        G4double innerRad = std::stod(tokens[1]);
        G4double outerRad = std::stod(tokens[2]);
        G4double height = std::stod(tokens[3]);
        G4double x = std::stod(tokens[4]);
        G4double y = std::stod(tokens[5]);
        G4double z = std::stod(tokens[6]);
        
        G4cout << "AddCargoCylinder: " << tokens[0] << " " << innerRad << " " << outerRad << " " << height 
               << " " << x << " " << y << " " << z << " " << tokens[7] << G4endl;
        
        AddCargoCylinder(tokens[0], 
                        innerRad*cm, outerRad*cm, height*cm,
                        x*cm, y*cm, z*cm,
                        tokens[7]);
    }
    catch (const std::exception& e) {
        G4cerr << "Error parsing parameters: " << e.what() << G4endl;
    }
}



G4VPhysicalVolume *RPCConstruction::Construct() {

    // Create world
    solidWorld = new G4Box("solidWorld", xWorldFull/2, yWorldFull/2, zWorldFull/2);
    logicWorld = new G4LogicalVolume(solidWorld, GetMaterialByName(worldMaterial), "logicWorld");
    G4UserLimits* worldLimit = new G4UserLimits();  
    worldLimit->SetMaxAllowedStep(surroundingStepLimit*mm);
    logicWorld->SetUserLimits(worldLimit);
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);  

    // Create cargo shapes
    for (auto &pair : cargoBoxes) {
        CargoBox &box = pair.second;
        G4Box *solidBox = new G4Box(box.name + "_solid", box.length/2, box.width/2, box.height/2);
        box.logicalVolume = new G4LogicalVolume(solidBox, box.material, box.name + "_logic");
        
        // Set step limits to ensure accurate simulation
        G4UserLimits* boxLimits = new G4UserLimits();  
        boxLimits->SetMaxAllowedStep(cargoStepLimit*mm);
        box.logicalVolume->SetUserLimits(boxLimits);
        
        new G4PVPlacement(0, box.position, box.logicalVolume, box.name + "_phys", logicWorld, false, 0, true);
    }
    
    for (auto &pair : cargoEllipsoids) {
        CargoEllipsoid &ellipsoid = pair.second;
        G4Ellipsoid *solidEllipsoid = new G4Ellipsoid(ellipsoid.name + "_solid", ellipsoid.xSemiAxis, ellipsoid.ySemiAxis, ellipsoid.zSemiAxis);
        ellipsoid.logicalVolume = new G4LogicalVolume(solidEllipsoid, ellipsoid.material, ellipsoid.name + "_logic");
        G4UserLimits* ellipsoidLimits = new G4UserLimits();
        ellipsoidLimits->SetMaxAllowedStep(cargoStepLimit*mm);
        ellipsoid.logicalVolume->SetUserLimits(ellipsoidLimits);
        new G4PVPlacement(0, ellipsoid.position, ellipsoid.logicalVolume, ellipsoid.name + "_phys", logicWorld, false, 0, true);
    }
    
    for (auto &pair : cargoCylinders) {
        CargoCylinder &cylinder = pair.second;
        G4Tubs *solidCylinder = new G4Tubs(cylinder.name + "_solid", cylinder.innerRadius, cylinder.outerRadius, cylinder.height / 2, 0, 360*deg);
        cylinder.logicalVolume = new G4LogicalVolume(solidCylinder, cylinder.material, cylinder.name + "_logic");
        G4UserLimits* cylinderLimits = new G4UserLimits();
        cylinderLimits->SetMaxAllowedStep(cargoStepLimit*mm);
        cylinder.logicalVolume->SetUserLimits(cylinderLimits);
        new G4PVPlacement(0, cylinder.position, cylinder.logicalVolume, cylinder.name + "_phys", logicWorld, false, 0, true);
    }

    // Create RPC "pixels"
    xPixelFull = pixelLength*cm;
    yPixelFull = pixelWidth*cm;
    zPixelFull = pixelThickness*mm;
    solidPixel = new G4Box("solidPixel", xPixelFull/2, yPixelFull/2, zPixelFull/2);
    logicPixel = new G4LogicalVolume(solidPixel, GetMaterialByName(rpcMaterial), "logicPixel");
    G4UserLimits* pixelLimit = new G4UserLimits();  
    pixelLimit->SetMaxAllowedStep(surroundingStepLimit*mm);
    logicPixel->SetUserLimits(pixelLimit);
    
    for (G4int i = 0; i < numY; i++) {
        for (G4int j = 0; j < numX; j++) {
            G4double xTrans = -(numX-1)/2.*xPixelFull + j * xPixelFull;
            G4double yTrans = -(numY-1)/2.*yPixelFull + i * yPixelFull;
            for (G4int k = 0; k < rpcHeights.size(); k++) {
                G4double zTrans = rpcHeights[k];
                G4Translate3D transRPC(G4ThreeVector(xTrans, yTrans, zTrans));
                physPixel = new G4PVPlacement(transRPC, logicPixel, "physPixel", logicWorld, false, numY*numX*k+numX*i+j, false);
            }
        }
    }

    return physWorld;
}

void RPCConstruction::ConstructSDandField() {

    EventAction *eventAction = static_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());

    SensitiveDetector *rpcdet = new SensitiveDetector("rpcdet", eventAction);
    logicPixel->SetSensitiveDetector(rpcdet);
    
    SensitiveDetector *worlddet = new SensitiveDetector("worlddet", eventAction);
    logicWorld->SetSensitiveDetector(worlddet);

    for (auto& pair : cargoBoxes) {
        SensitiveDetector *boxdet = new SensitiveDetector(pair.first, eventAction);
        pair.second.logicalVolume->SetSensitiveDetector(boxdet);
    }
    for (auto& pair : cargoEllipsoids) {
        SensitiveDetector *boxdet = new SensitiveDetector(pair.first, eventAction);
        pair.second.logicalVolume->SetSensitiveDetector(boxdet);
    }
    for (auto& pair : cargoCylinders) {
        SensitiveDetector *boxdet = new SensitiveDetector(pair.first, eventAction);
        pair.second.logicalVolume->SetSensitiveDetector(boxdet);
    }
}


