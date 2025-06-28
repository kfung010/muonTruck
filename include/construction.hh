#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
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
#include <map>

struct CargoBox {
    G4String name;
    G4double width;
    G4double length;
    G4double height;
    G4ThreeVector position;
    G4Material* material;
    G4LogicalVolume* logicalVolume = nullptr;
};

struct CargoEllipsoid {
    G4String name;
    G4double xSemiAxis;    // Half-length along X
    G4double ySemiAxis;    // Half-length along Y
    G4double zSemiAxis;    // Half-length along Z
    G4ThreeVector position;
    G4Material* material;
    G4LogicalVolume* logicalVolume = nullptr;
};

struct CargoCylinder {
    G4String name;
    G4double innerRadius;  // Set to 0 for solid cylinder
    G4double outerRadius;
    G4double height;       // Full height (along Z-axis)
    G4ThreeVector position;
    G4Material* material;
    G4LogicalVolume* logicalVolume = nullptr;
};



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
        
        G4double getzHighestRPC() const { return *std::max_element(rpcHeights.begin(), rpcHeights.end()) + zPixelFull/2; }
        
        void ClearRPCHeights();
        void AddRPCHeight(G4double height);
        
        void AddCargoBox(const G4String& params);
        void AddCargoEllipsoid(const G4String& params);
        void AddCargoCylinder(const G4String& params);

        void SetupCargoLayout();
        
        void ClearCargo();

    private:

        G4double xWorldFull = 50*m; 
        G4double yWorldFull = 50*m; 
        G4double zWorldFull = 50*m; 

        G4double xPixelFull, yPixelFull, zPixelFull;
        std::vector<G4double> rpcHeights = {-290*cm, -270*cm, -250*cm, -230*cm, 230*cm, 250*cm, 270*cm, 290*cm};
        G4int numX = 1;
        G4int numY = 1;

        G4Box *solidWorld, *solidPixel;
        G4LogicalVolume *logicWorld, *logicPixel;
        G4VPhysicalVolume *physWorld, *physPixel;
        
        G4double pixelThickness = 1;
        G4double pixelLength = 1000;
        G4double pixelWidth = 1000;
        G4double cargoStepLimit = 0.01;
        G4double surroundingStepLimit = 100000;
        

        void DefineMaterials();
        G4String worldMaterial = "air";
        G4String rpcMaterial = "air";
        G4Material *worldMat, *rpcMat;
        
        G4Material *ammoniumNitrate, *H2O, *air, *vacuum, *C2H2F4, *iC4H10, *RPCgas;
        G4Material *lithium, *aluminium, *copper, *lead, *tungsten, *carbon;

        virtual void ConstructSDandField();
        
        G4GenericMessenger *fMessenger;
        
        std::map<G4String, CargoBox> cargoBoxes;
        std::map<G4String, CargoEllipsoid> cargoEllipsoids;
        std::map<G4String, CargoCylinder> cargoCylinders;
        
        G4Material* GetMaterialByName(const G4String& materialName) const {
            if (materialName == "ammoniumNitrate") return ammoniumNitrate;
            if (materialName == "water") return H2O;
            if (materialName == "air") return air;
            if (materialName == "vacuum") return vacuum;
            if (materialName == "C2H2F4") return C2H2F4;
            if (materialName == "iC4H10") return iC4H10;
            if (materialName == "RPCgas") return RPCgas;
            if (materialName == "lithium") return lithium;
            if (materialName == "aluminium") return aluminium;
            if (materialName == "copper") return copper;
            if (materialName == "lead") return lead;
            if (materialName == "tungsten") return tungsten;
            if (materialName == "carbon") return carbon;
            throw std::runtime_error("Invalid material: " + materialName + ". Program ended with an error.");
        }
        
        void AddCargoBox(const G4String& name, G4double width, G4double height, G4double length, 
                         G4double posX, G4double posY, G4double posZ, const G4String& materialName);
        
        void AddCargoEllipsoid(const G4String& name, 
                       G4double xSemiAxis, G4double ySemiAxis, G4double zSemiAxis,
                       G4double posX, G4double posY, G4double posZ, const G4String& materialName);

        void AddCargoCylinder(const G4String& name, 
                              G4double innerRadius, G4double outerRadius, G4double height,
                              G4double posX, G4double posY, G4double posZ, const G4String& materialName);
};

#endif
