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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(G4String FilenameRoot)
 : G4UserRunAction(),
   NameRootdat(FilenameRoot)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  //analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
  //analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
  //analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  //analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

  // Creating ntuple
  //
  // Creating ntuple vacstep1 , id=0
  //
  analysisManager->CreateNtuple("Spectrum", "vacstep3");
  analysisManager->CreateNtupleIColumn("pdg");
  analysisManager->CreateNtupleDColumn("E");
  /*analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("startx");
  analysisManager->CreateNtupleDColumn("starty");
  analysisManager->CreateNtupleDColumn("startz");
  analysisManager->CreateNtupleDColumn("px");
  analysisManager->CreateNtupleDColumn("py");
  analysisManager->CreateNtupleDColumn("pz");
  analysisManager->CreateNtupleDColumn("Polx");
  analysisManager->CreateNtupleDColumn("Poly");
  analysisManager->CreateNtupleDColumn("Polz");
  analysisManager->CreateNtupleDColumn("TrackID");
  analysisManager->CreateNtupleDColumn("ParentID");
  analysisManager->CreateNtupleDColumn("EventID");*/
  analysisManager->FinishNtuple();
/*
  // Creating ntuple vacstep2 , id=1
  //
  analysisManager->CreateNtuple("bremssim2", "vacstep2");
  analysisManager->CreateNtupleIColumn("pdg");
  analysisManager->CreateNtupleDColumn("E");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("startx");
  analysisManager->CreateNtupleDColumn("starty");
  analysisManager->CreateNtupleDColumn("startz");
  analysisManager->CreateNtupleDColumn("px");
  analysisManager->CreateNtupleDColumn("py");
  analysisManager->CreateNtupleDColumn("pz");
  analysisManager->CreateNtupleDColumn("Polx");
  analysisManager->CreateNtupleDColumn("Poly");
  analysisManager->CreateNtupleDColumn("Polz");
  analysisManager->CreateNtupleDColumn("TrackID");
  analysisManager->CreateNtupleDColumn("ParentID");
  analysisManager->CreateNtupleDColumn("EventID");
  analysisManager->FinishNtuple();

   // Creating ntuple Detector , id=2
  //
  analysisManager->CreateNtuple("bremssim3", "Detector");
  analysisManager->CreateNtupleIColumn("pdg");
  analysisManager->CreateNtupleDColumn("E");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("startx");
  analysisManager->CreateNtupleDColumn("starty");
  analysisManager->CreateNtupleDColumn("startz");
  analysisManager->CreateNtupleDColumn("px");
  analysisManager->CreateNtupleDColumn("py");
  analysisManager->CreateNtupleDColumn("pz");
  analysisManager->CreateNtupleDColumn("Polx");
  analysisManager->CreateNtupleDColumn("Poly");
  analysisManager->CreateNtupleDColumn("Polz");
  analysisManager->CreateNtupleDColumn("TrackID");
  analysisManager->CreateNtupleDColumn("ParentID");
  analysisManager->CreateNtupleDColumn("EventID");
  analysisManager->FinishNtuple();
*/
  // Creating ntuple Detector Energie, id=3
  analysisManager->CreateNtuple("B4", "EdepandTrackL");
  analysisManager->CreateNtupleDColumn("Eabs");
  //analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Labs");
  //analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
 // G4String fileName = "results";
  analysisManager->OpenFile(NameRootdat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  /*
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }

    G4cout << " EAbs : mean = "
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

    G4cout << " EGap : mean = "
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

    G4cout << " LAbs : mean = "
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

    G4cout << " LGap : mean = "
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
  }
*/
  // save histograms & ntuple
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
