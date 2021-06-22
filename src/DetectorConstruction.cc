//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4PolarizationManager.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fMagnetPV(nullptr),
   fDetectorLV(nullptr),
   fdetmat(0),
   fCheckOverlaps(true)

{ //default parameters
	magthick =  15.*cm;
    absthick = 2.*mm;
  //materials
   DefineMaterials();
   SetDetectorMaterial("LeadGlass");

	fMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
	delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  G4int ncomponents, natoms;
  G4String name, symbol;
  G4double fractionmass;

  //
  // define Elements
  //


  G4Element* I       = new G4Element(name="Iodine"  ,symbol="I"   , z= 53., a = 126.904*g/mole);
  G4Element* Cs      = new G4Element(name="Cesium"  ,symbol="Cs"  , z= 55., a = 132.905*g/mole);
  G4Element* O       = new G4Element(name="Oxygen"  ,symbol="O"   , z= 8  , a =  16.000*g/mole);
  G4Element* As      = new G4Element(name="Arsenic" ,symbol="As"  , z= 33 , a =   74.922*g/mole);


  G4Element* C      = new G4Element(name="Carbon" ,symbol="C"  , z= 6 , a =   12*g/mole);
  G4Element* H     = new G4Element(name="Hydrogen" ,symbol="H"  , z= 1 , a =   2*g/mole);


  // definde materials
  density = 3.74*g/cm3;
  G4Material* As2O3 = new G4Material(name="Arsenictrioxide", density, ncomponents=2);
  As2O3->AddElement(As, natoms=2);
  As2O3->AddElement(O , natoms=3);

  //material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_W");
  G4Material* PbO = nistManager->FindOrBuildMaterial("G4_LEAD_OXIDE");
  G4Material* SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* K2O = nistManager->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
  // CsI
  density = 4.53*g/cm3;
  G4Material* CsI = new G4Material(name="CesiumIodide", density, ncomponents=2);
  CsI->AddElement(Cs, natoms=5);
  CsI->AddElement(I, natoms=5);

  // Lead Glass TF1-000
  G4Material* LeadGlass = new G4Material("LeadGlass", density= 3.860*g/cm3, ncomponents=4);
  LeadGlass->AddMaterial(PbO   , fractionmass=0.512);
  LeadGlass->AddMaterial(SiO2  , fractionmass=0.413);
  LeadGlass->AddMaterial(K2O   , fractionmass=0.07);
  LeadGlass->AddMaterial(As2O3 , fractionmass=0.005);
  //*************************************************************************************
  // This part below is needed to switsh on the scintillation and optical photons process
  //                      (not relevant for the asymetry)
  //*************************************************************************************
  const G4int NUMENTRIES = 2;
  G4double CsI_PP[NUMENTRIES]    = { 171.2068*eV, 331*eV };
  G4double CsI_SCINT[NUMENTRIES] = { 1.0, 1.0 };
  G4double CsI_RIND[NUMENTRIES]  = { 1.57, 1.57 };
  G4double CsI_ABSL[NUMENTRIES]  = { 1500*cm, 2500*cm};
  G4MaterialPropertiesTable *CsI_mt = new G4MaterialPropertiesTable();
  CsI_mt->AddProperty("SCINTILLATION", CsI_PP, CsI_SCINT, NUMENTRIES);
  CsI_mt->AddProperty("RINDEX",        CsI_PP, CsI_RIND,  NUMENTRIES);
  CsI_mt->AddProperty("ABSLENGTH",     CsI_PP, CsI_ABSL,  NUMENTRIES);
  CsI->SetMaterialPropertiesTable(CsI_mt);


// this is are the optical properties of quartz just to test lead Glas (taken from Quartz Cherenkov Simulation)
  const G4int nSpectRI=25;
  const G4int nSpectAbs=9;
  G4double EnergySpectrosil[nSpectRI]={1.38 *eV, 1.455*eV, 1.755*eV, 1.889*eV,
                                       1.926*eV, 1.959*eV, 2.104*eV, 2.110*eV,
                                       2.270*eV, 2.550*eV, 2.845*eV, 3.064*eV,
                                       3.397*eV, 3.709*eV, 3.965*eV,  4.886*eV,
                                       4.993*eV, 4.999*eV, 5.417*eV,  5.780*eV,
                                       6.011*eV, 6.383*eV, 6.411*eV, 6.424*eV, 6.7*eV};
  G4double RISpectrosil[nSpectRI]={1.45181, 1.45247, 1.45515, 1.45637, 1.4567 ,
                                   1.45702, 1.4584 , 1.45846, 1.46008, 1.46313,
                                   1.46669, 1.46962, 1.47454, 1.47975, 1.48447,
                                   1.50547, 1.50838, 1.50855, 1.52109, 1.53365,
                                   1.54259, 1.55884, 1.56014, 1.56077, 1.57495};
  G4double ESpectrosil[nSpectAbs] ={1.38*eV , 3.1*eV ,  3.35*eV , 4.0*eV , 4.43*eV , 4.96*eV , 5.64*eV , 6.53*eV , 6.7*eV};
  G4double ABSpectrosil[nSpectAbs]={134.8*mm, 131.2*mm, 130.1*mm, 124.8*mm, 119.9*mm, 112.6*mm, 106.0*mm, 94.9*mm, 1.6*mm};
  G4MaterialPropertiesTable* MPT_Quartz = new G4MaterialPropertiesTable();


  MPT_Quartz->AddProperty ("RINDEX",  EnergySpectrosil, RISpectrosil, nSpectRI);
  MPT_Quartz->AddProperty ("ABSLENGTH", ESpectrosil, ABSpectrosil, nSpectAbs);

  LeadGlass->SetMaterialPropertiesTable (MPT_Quartz);





  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double worldSizeXY = 3.*m;
  G4double worldSizeZ  = 4.*m;
  //G4double absthick = 2.*mm;
  G4double absrad=25.*mm;
  //G4double magthick =  15.*cm;
  G4double vacthick  = 1.*mm;
  G4double gap1=12.5*mm;
  G4double gap2=10.*mm;
  G4double gap3=10.*mm; //gap for detector
  G4double detthick = 45.*cm; // detector size
  G4double detx = 3.8 *cm;
  G4double dety = 3.8 *cm;
  //G4double vac3x = 40. *cm;
  //G4double vac3y = 40. *cm;

  G4double vac3ir = 2.8 *cm;
  G4double vac3ar = 2.9*cm;




  auto ZposAbs=0*cm;
  auto ZposMag=absthick/2.0 + gap1 + magthick/2.0;
  auto ZposVac1=absthick/2.0 +gap1-1.*mm-vacthick/2.0;
  auto ZposVac2=absthick/2.0+gap1+magthick+gap2+vacthick/2.0;
  auto ZposDet=absthick/2.0+gap1+magthick+gap2+vacthick+gap3+detthick/2.0; // z position of the detektor
  //auto ZposVac3=absthick/2.0+gap1+magthick+gap2+vacthick/2.0+detthick+gap2+vacthick/2.0;

  // Get materials
  auto worldMat = G4Material::GetMaterial("Galactic");
  //auto absMat = G4Material::GetMaterial("G4_W");
  auto absMat = G4Material::GetMaterial("Galactic");
//  auto magMat = G4Material::GetMaterial("G4_Fe");
  auto magMat = G4Material::GetMaterial("Galactic");
  if ( ! worldMat || ! absMat || ! magMat ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 worldMat,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Absorber
  //
  auto absorberS
    = new G4Tubs("Absorber", //Name
                0.,         // inner radius
                absrad,     // outer radius
                absthick/2., // half length in z
                0.0*rad,    // starting phi angle
                twopi*rad); // angle of the segment

  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,     // its solid
                 absMat,  // its material
                 "Absorber");   // its name
  fAbsorberPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 absorberLV,          // its logical volume
                 "Absorber",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


   //
   //vacuum step 1
   //
   auto
   fVacStepS1 = new G4Tubs("VacStep",  //Name
                               0.,         // inner radius
                               absrad,     // outer radius
                               vacthick/2., // half length in z
                               0.0*deg,    // starting phi angle
                               360.0*deg); // angle of the segment

  auto
   fVacStepLV1 = new G4LogicalVolume(fVacStepS1,    //its solid
                                        worldMat,    //its material
                                        "VacStep");  //its name

   fVacStepPV1 = new G4PVPlacement(0,                   //no rotation
                        G4ThreeVector(0.,0., ZposVac1),    //its position
                                fVacStepLV1,            //its logical volume
                                "VacStep",                 //its name
                                worldLV,               //its mother
                                false,                     //no boolean operat
                                0);                        //copy number

  //
  // Magnet
  //
 auto magnetS
   = new G4Tubs("magnetS", //Name
               0.,         // inner radius
               absrad,     // outer radius
               magthick/2., // half length in z
               0.0*rad,    // starting phi angle
               twopi*rad); // angle of the segment

  auto magnetLV
    = new G4LogicalVolume(
                 magnetS,        // its solid
                 magMat, // its material
                 "magnetLV");          // its name

  fMagnetPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., ZposMag), // its position
                 magnetLV,       // its logical volume
                 "magnetPV",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

 auto polvec=G4ThreeVector(0.,0.,1.0);
 G4PolarizationManager * polMgr = G4PolarizationManager::GetInstance();
 polMgr->SetVolumePolarization(magnetLV, polvec);

 //
 //vacuum step 2
 //
auto fVacStepS2 = new G4Tubs("VacStep2",  //Name
                             0.,         // inner radius
                             absrad,     // outer radius
                             vacthick/2., // half length in z
                             0.0*deg,    // starting phi angle
                             360.0*deg); // angle of the segment

 auto fVacStepLV2 = new G4LogicalVolume(fVacStepS2,    //its solid
                                      worldMat,    //its material
                                      "VacStep2");       //its name

 fVacStepPV2 = new G4PVPlacement(0,                   //no rotation
                      G4ThreeVector(0.,0., ZposVac2),    //its position
                              fVacStepLV2,            //its logical volume
                              "VacStep2",                 //its name
                              worldLV,               //its mother
                              false,                     //no boolean operat
                              0);                        //copy number

 //
 // Detector
 //
 //Detector Box shape as in E166
auto fDetectorS= new G4Box("Detector",  //Name
                             detx/2.,   // x size
                             dety/2.,     // y size
                             detthick/2.); // z size

/*
auto fDetectorS= new G4Tubs("Detector",  //Name
                             0.,         // inner radius
                             absrad,     // outer radius
                             detthick/2., // half length in z
                             0.0*deg,    // starting phi angle
                             360.0*deg); // angle of the segment
*/
 fDetectorLV = new G4LogicalVolume(fDetectorS,    //its solid
                                      fdetmat,    //its material
                                      "Detector");       //its name

 fDetectorPV = new G4PVPlacement(0,                   //no rotation
                      G4ThreeVector(0.,0., ZposDet),    //its position
                              fDetectorLV,            //its logical volume
                              "Detector",                 //its name
                              worldLV,               //its mother
                              false,                     //no boolean operat
                              0);                        //copy number

//
 //vacuum step 3
 /*
auto fVacStepS3 = new G4Box("VacStep3",  //Name
                             vac3x/2.,
                             vac3y/2,
                             vacthick/2.);

 auto fVacStepLV3 = new G4LogicalVolume(fVacStepS3,    //its solid
                                      worldMat,    //its material
                                      "VacStep3");       //its name

 fVacStepPV3 = new G4PVPlacement(0,                   //no rotation
                      G4ThreeVector(0.,0., ZposVac3),    //its position
                              fVacStepLV3,            //its logical volume
                              "VacStep3",                 //its name
                              worldLV,               //its mother
                              false,                     //no boolean operat
                              0);                        //copy number

*/

//
//vacuum step 3 tube
//
auto fVacStepS3 = new G4Tubs("VacStep3",  //Name
                            vac3ir,         // inner radius
                            vac3ar,     // outer radius
                            detthick/2., // half length in z
                            0.0*deg,    // starting phi angle
                            360.0*deg); // angle of the segment

auto fVacStepLV3 = new G4LogicalVolume(fVacStepS3,    //its solid
                                     worldMat,    //its material
                                     "VacStep3");       //its name

fVacStepPV3 = new G4PVPlacement(0,                   //no rotation
                     G4ThreeVector(0.,0.,ZposDet),    //its position
                             fVacStepLV3,            //its logical volume
                             "VacStep3",                 //its name
                             worldLV,               //its mother
                             false,                     //no boolean operat
                             0);                        //copy number






  //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "\n The  WORLD "<< G4endl
    << "\t material: "<< worldMat->GetName()<< G4endl
    << "\t Size in Z: "<< G4BestUnit(worldSizeZ,"Length")<< G4endl
    << "\t Size in XY: "<< G4BestUnit(worldSizeXY,"Length")<< G4endl

    << "\n The Reconversion Target "<< G4endl
    << "\t material: "<< absMat->GetName()<< G4endl
    << "\t thickness: "<< G4BestUnit(absthick,"Length")<< G4endl
    << "\t radius: "<< G4BestUnit(absrad,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposAbs,"Length")<< G4endl

    << "\n The First Vacuum Step "<< G4endl
    << "\t thickness: "<< G4BestUnit(vacthick,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposVac1,"Length")<< G4endl

    << "\n The Magnetized Iron Block "<< G4endl
    << "\t material: "<< magMat->GetName()<< G4endl
    << "\t thickness: "<< G4BestUnit(magthick,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposMag,"Length")<< G4endl
    //<< "\t Z-component of B-field: "<< G4BestUnit(magvec[2],"Magnetic flux density")<< G4endl
    << "\t Z-component of Polarization: "<< polvec[2] << G4endl

    << "\n The Second Vacuum Step "<< G4endl
    << "\t thickness: "<< G4BestUnit(vacthick,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposVac2,"Length")<< G4endl

    << "\n The Detector "<< G4endl
    << "\t material: "<< fdetmat->GetName()<< G4endl
    << "\t thickness: "<< G4BestUnit(detthick,"Length")<< G4endl
    << "\t radius: "<< G4BestUnit(absrad,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposDet,"Length")<< G4endl

    << "\n The Fourth Vacuum Step "<< G4endl
    << "\t thickness: "<< G4BestUnit(detthick,"Length")<< G4endl
    << "\t Inner radius: "<< G4BestUnit(vac3ir,"Length")<< G4endl
    << "\t Outer radius: "<< G4BestUnit(vac3ar,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposDet,"Length")<< G4endl
/*
    << "\n The Third  Vacuum Step "<< G4endl
    << "\t thickness: "<< G4BestUnit(vacthick,"Length")<< G4endl
    << "\t Z-position of centre: "<< G4BestUnit(ZposVac3,"Length")<< G4endl
*/

    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto MagVisAtt= new G4VisAttributes(G4Colour(1.0,0.5,1.0));
  MagVisAtt->SetVisibility(true);
  magnetLV->SetVisAttributes(MagVisAtt);

  auto AbsVisAtt= new G4VisAttributes(G4Colour(0.5,1.0,1.0));
  AbsVisAtt->SetVisibility(true);
  absorberLV->SetVisAttributes(AbsVisAtt);

  auto DetVisAtt= new G4VisAttributes(G4Colour(0.5,1.0,0.5));
  DetVisAtt->SetVisibility(true);
  fDetectorLV->SetVisAttributes(DetVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}


void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fdetmat != pttoMaterial) {
    fdetmat = pttoMaterial;
    if(fDetectorLV) fDetectorLV->SetMaterial(fdetmat);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void DetectorConstruction::SetTargetThicknes(G4double val)
{
  // change Target thickness (Tungstenabsorber)
  absthick = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Thicknes of the Thunstentarget :"<< G4BestUnit(GetTargetTicknes(),"Length") << G4endl;
}

void DetectorConstruction::SetAbsorberThicknes(G4double val)
{
  // change Absorber thickness (Ironabsorber)
  magthick = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Thicknes of the Ironabsorber :" << G4BestUnit(GetAbsorberTicknes(),"Length")  << G4endl;
}

 G4VPhysicalVolume* DetectorConstruction::UpdateGeometry()
{
	return DefineVolumes();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
