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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Analysis.hh"
// all header files below are coppied from OpNoviceStackingAction.cc
#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
// heder files below are coppied from QuaSiPMT.cc
#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  //auto particle = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  auto aTrack = step->GetTrack();

  if ( volume == fDetConstruction->GetVacStep3PV() /*&& particle == 22*/) {
       // get analysis manager
     // if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {

       auto analysisManager = G4AnalysisManager::Instance();

       // fill ntuple id=0
       analysisManager->FillNtupleIColumn(0,0, step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());

       analysisManager->FillNtupleDColumn(0,1,step->GetPostStepPoint()->GetTotalEnergy()/eV );
/*
       analysisManager->FillNtupleDColumn(0,2, step->GetPostStepPoint()->GetPosition().x());
       analysisManager->FillNtupleDColumn(0,3, step->GetPostStepPoint()->GetPosition().y());
       analysisManager->FillNtupleDColumn(0,4, step->GetPostStepPoint()->GetPosition().z());

       analysisManager->FillNtupleDColumn(0,5, step->GetTrack()->GetVertexPosition().x());
       analysisManager->FillNtupleDColumn(0,6, step->GetTrack()->GetVertexPosition().y());
       analysisManager->FillNtupleDColumn(0,7, step->GetTrack()->GetVertexPosition().z());

       analysisManager->FillNtupleDColumn(0,8, step->GetPostStepPoint()->GetMomentumDirection().x());
       analysisManager->FillNtupleDColumn(0,9, step->GetPostStepPoint()->GetMomentumDirection().y());
       analysisManager->FillNtupleDColumn(0,10, step->GetPostStepPoint()->GetMomentumDirection().z());

       analysisManager->FillNtupleDColumn(0,11, step->GetTrack()->GetPolarization().x());
       analysisManager->FillNtupleDColumn(0,12, step->GetTrack()->GetPolarization().y());
       analysisManager->FillNtupleDColumn(0,13, step->GetTrack()->GetPolarization().z());

       analysisManager->FillNtupleDColumn(0,14, step->GetTrack()->GetTrackID());
       analysisManager->FillNtupleDColumn(0,15, step->GetTrack()->GetParentID());
       analysisManager->FillNtupleDColumn(0,16, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
*/
       analysisManager->AddNtupleRow(0);

  //G4cout<< " This part of the code you are currently testing is executed"  << G4endl;
     }
   // }
/*
     if ( volume == fDetConstruction->GetVacStep2PV() ) {

          auto analysisManager = G4AnalysisManager::Instance();
         //~//~//~//
         analysisManager->FillNtupleIColumn(1,0, step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());

         analysisManager->FillNtupleDColumn(1,1, step->GetPostStepPoint()->GetKineticEnergy()/MeV);

         analysisManager->FillNtupleDColumn(1,2, step->GetPostStepPoint()->GetPosition().x());
         analysisManager->FillNtupleDColumn(1,3, step->GetPostStepPoint()->GetPosition().y());
         analysisManager->FillNtupleDColumn(1,4, step->GetPostStepPoint()->GetPosition().z());

         analysisManager->FillNtupleDColumn(1,5, step->GetTrack()->GetVertexPosition().x());
         analysisManager->FillNtupleDColumn(1,6, step->GetTrack()->GetVertexPosition().y());
         analysisManager->FillNtupleDColumn(1,7, step->GetTrack()->GetVertexPosition().z());

         analysisManager->FillNtupleDColumn(1,8, step->GetPostStepPoint()->GetMomentumDirection().x());
         analysisManager->FillNtupleDColumn(1,9, step->GetPostStepPoint()->GetMomentumDirection().y());
         analysisManager->FillNtupleDColumn(1,10, step->GetPostStepPoint()->GetMomentumDirection().z());

         analysisManager->FillNtupleDColumn(1,11, step->GetTrack()->GetPolarization().x());
         analysisManager->FillNtupleDColumn(1,12, step->GetTrack()->GetPolarization().y());
         analysisManager->FillNtupleDColumn(1,13, step->GetTrack()->GetPolarization().z());

         analysisManager->FillNtupleDColumn(1,14, step->GetTrack()->GetTrackID());
         analysisManager->FillNtupleDColumn(1,15, step->GetTrack()->GetParentID());
         analysisManager->FillNtupleDColumn(1,16, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());

         analysisManager->AddNtupleRow(1);
	 }

     if ( volume == fDetConstruction->GetDetectorPV() ) {

          auto analysisManager = G4AnalysisManager::Instance();
         //~//~//~//
         analysisManager->FillNtupleIColumn(2,0, step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());

         analysisManager->FillNtupleDColumn(2,1, step->GetPostStepPoint()->GetKineticEnergy()/MeV);

         analysisManager->FillNtupleDColumn(2,2, step->GetPostStepPoint()->GetPosition().x());
         analysisManager->FillNtupleDColumn(2,3, step->GetPostStepPoint()->GetPosition().y());
         analysisManager->FillNtupleDColumn(2,4, step->GetPostStepPoint()->GetPosition().z());

         analysisManager->FillNtupleDColumn(2,5, step->GetTrack()->GetVertexPosition().x());
         analysisManager->FillNtupleDColumn(2,6, step->GetTrack()->GetVertexPosition().y());
         analysisManager->FillNtupleDColumn(2,7, step->GetTrack()->GetVertexPosition().z());

         analysisManager->FillNtupleDColumn(2,8, step->GetPostStepPoint()->GetMomentumDirection().x());
         analysisManager->FillNtupleDColumn(2,9, step->GetPostStepPoint()->GetMomentumDirection().y());
         analysisManager->FillNtupleDColumn(2,10, step->GetPostStepPoint()->GetMomentumDirection().z());

         analysisManager->FillNtupleDColumn(2,11, step->GetTrack()->GetPolarization().x());
         analysisManager->FillNtupleDColumn(2,12, step->GetTrack()->GetPolarization().y());
         analysisManager->FillNtupleDColumn(2,13, step->GetTrack()->GetPolarization().z());

         analysisManager->FillNtupleDColumn(2,14, step->GetTrack()->GetTrackID());
         analysisManager->FillNtupleDColumn(2,15, step->GetTrack()->GetParentID());
         analysisManager->FillNtupleDColumn(2,16, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());

         analysisManager->AddNtupleRow(2);
    //G4cout<< " This part of the code you are currently testing is executed"  << G4endl;
       }
*/
// Collect energy and track length step by step



  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  //if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
   // stepLength = step->GetStepLength();
  //h}
  stepLength = step->GetStepLength();

  if ( volume == fDetConstruction->GetDetectorPV()  ) {
    fEventAction->AddAbs(edep,stepLength);
  }
  /*
  if ( volume == fDetConstruction->GetGapPV() ) {
    fEventAction->AddGap(edep,stepLength);
  }*/
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
