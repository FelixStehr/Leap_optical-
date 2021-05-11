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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fDirectory = new G4UIdirectory("/leapsim/");
  fDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/leapsim/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fDetMatCmd = new G4UIcmdWithAString("/leapsim/det/setDetMat",this);
  fDetMatCmd->SetGuidance("Select Material of the Detector.");
  fDetMatCmd->SetParameterName("choice",false);
  fDetMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
/*
  fChamMatCmd = new G4UIcmdWithAString("/B2/det/setChamberMaterial",this);
  fChamMatCmd->SetGuidance("Select Material of the Chamber.");
  fChamMatCmd->SetParameterName("choice",false);
  fChamMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/B2/det/stepMax",this);
  fStepMaxCmd->SetGuidance("Define a step max");
  fStepMaxCmd->SetParameterName("stepMax",false);
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_Idle);
  */
  
  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/leapsim/det/setAbsThick",this);
  fAbsThickCmd->SetGuidance("Define Absorber Thicknes");
  fAbsThickCmd->SetParameterName("AbsThick",false);
  fAbsThickCmd->SetUnitCategory("Length");
  fAbsThickCmd->AvailableForStates(G4State_Idle);
  
  fTargThickCmd= new G4UIcmdWithADoubleAndUnit("/leapsim/det/setTargThick",this);
  fTargThickCmd->SetGuidance("Define Target Thicknes");
  fTargThickCmd->SetParameterName("TargThick",false);
  fTargThickCmd->SetUnitCategory("Length");
  fTargThickCmd->AvailableForStates(G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/leapsim/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetMatCmd;
  //delete fChamMatCmd;
  //delete fStepMaxCmd;
  //delete fB2Directory;
  delete fTargThickCmd;
  delete fAbsThickCmd;
  delete fDetDirectory;
  delete fDirectory;
  delete UpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if ( command == fDetMatCmd )
   {fDetectorConstruction->SetDetectorMaterial(newValue);}
   
  if( command == fTargThickCmd )
   { fDetectorConstruction->SetTargetThicknes(fTargThickCmd->GetNewDoubleValue(newValue));}

  if( command == fAbsThickCmd )
   { fDetectorConstruction->SetAbsorberThicknes(fAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { fDetectorConstruction->UpdateGeometry(); }
   
  //if( command == fStepMaxCmd ) {
   // fDetectorConstruction
   //   ->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
