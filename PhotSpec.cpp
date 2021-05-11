#include<algorithm>
#include<iostream>
#include<vector>
#include "TTree.h"
#include "TFile.h"
using namespace std;




void PhotSpec(){
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int add_plot=0;

// get pointers to the trees
    // get file of magnetic field parallel to polarization and propagation direction of incoming e-
    TFile *fBz = new TFile("results.root");
    //TFile *f2 = new TFile("40MeV.root");
    // pointer to ttree of sensor after absorber
    //TTree *TBz1=(TTree*)fBz->Get("bremssim1"); //1
    //pointer to ttree of sensor after magnet
    TTree *TBz2=(TTree*)fBz->Get("Spectrum"); //2
    //pointer to ttree of sensor after magnet
    //TTree *TBz3=(TTree*)fBz->Get("bremssim3"); //3
    //pointer to ttree of deposited energy
   // TTree *TBz4=(TTree*)fBz->Get("B4"); //3
   // TTree *T24=(TTree*)f2->Get("B4"); //3
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//    cout<<PDG2[500]<<"\t"<<E2[500]<<"\t"<<TID2[500]<<"\t"<<PID2[500]<<"\t"<<EID2[500]<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Array with bins for energy binning
int nbins = 40;
double startval = 0.1;
double endval =60;
double ebins[nbins+1];
for (int i=0; i<=nbins; i++)
{
  ebins[i]=exp(log(startval)+(log(endval)-log(startval))/nbins*double(i));
  //cout<<ebins[i]<<endl;
}

int bins = 20;
double start = 0;
double end =60000000;
//double end = 7

TH1::SetDefaultSumw2();

//Create a new Canvas at random Position
TCanvas *C1 = new TCanvas("C1", // canvas name
                          "C1", // canvas title
                          1200, // canvas size in pixels along x
                          600); // canvas size in pixels along y

//Creates a 1D-histogram with floating numbers
  TH1F *PhoSpec= new TH1F("PhoSpec", //Name
                          "E_{Beam} = 60 MeV;Photon Energy / eV; Counts ",// title; axes label
                        //  nbins, ebins);// xbins
                         bins , start , end);
  TBz2->Draw("E>>PhoSpec");
  //TBz2->Draw("E>>PhoSpec");
  PhoSpec->SetLineWidth(2);
  //PhoSpec->SetMarkerStyle(kFullCircle);
  PhoSpec->GetXaxis()->SetTitleSize(0.05);
  PhoSpec->GetYaxis()->SetTitleSize(0.05);
  PhoSpec->GetXaxis()->SetLabelSize(0.05);
  PhoSpec->GetYaxis()->SetLabelSize(0.05);
  PhoSpec->GetYaxis()->SetTitleOffset(1.0);
gStyle->SetOptStat(1110);
C1->Print("PhoSpec.png");



/*


// First create histogram of Photons after absorber

//Create a new Canvas at random Position
TCanvas *C1 = new TCanvas("C1", // canvas name
                          "C1", // canvas title
                          900, // canvas size in pixels along x
                          600); // canvas size in pixels along y

//Creates a 1D-histogram with floating numbers
  TH1F *InP = new TH1F("InP", //Name
                          "100k incident electrons, S_{in}=(0,0,1), parallel magnetic field;Kinetic Energy / MeV; Number of Incident Photons",// title; axes label
                          nbins, // number of bins X
                          ebins);// xbins

  TBz1->Draw("E>>InP", // varexp; "e1" produces TH1F, "e1:e2" unbinned 2D scatter plot
          "pdg==22 & z<13.0 ", // selection; if boolean expression is true, hist is filled with a weight = value
          "E");
gStyle->SetOptStat(11);
C1->Print("gamma_InP.png");


// Create histogram of Photons after magnet

//Create a new Canvas at random Position
TCanvas *C2 = new TCanvas("C2", // canvas name
                          "C2", // canvas title
                          900, // canvas size in pixels along x
                          600); // canvas size in pixels along y

//Creates a 1D-histogram with floating numbers
  TH1F *TransP = new TH1F("TransP", //Name
                          "100k incident electrons, S_{in}=(0,0,1), parallel magnetic field;Kinetic Energy / MeV; Number of Incident Photons",// title; axes label
                          nbins, // number of bins X
                          ebins);// xbins

  TBz2->Draw("E>>TransP");

gStyle->SetOptStat(11);
C2->Print("gamma_TransP.png");

// Create histogram of Photons in the Detector

//Create a new Canvas at random Position
TCanvas *C3 = new TCanvas("C3", // canvas name
                          "C3", // canvas title
                          900, // canvas size in pixels along x
                          600); // canvas size in pixels along y

//Creates a 1D-histogram with floating numbers
  TH1F *DetP = new TH1F("DetP", //Name
                          "100k incident electrons, S_{in}=(0,0,1), parallel magnetic field;Kinetic Energy / MeV; Number of Incident Photons",// title; axes label
                          nbins, // number of bins X
                          ebins);// xbins

  TBz3->Draw("E>>DetP");

gStyle->SetOptStat(11);
C3->Print("gamma_DetP.png");

// Create histogram of Photons in the Detector2

//Create a new Canvas at random Position
TCanvas *C4 = new TCanvas("C4", // canvas name
                          "C4", // canvas title
                          900, // canvas size in pixels along x
                          600); // canvas size in pixels along y

  TBz3->Draw("E");

gStyle->SetOptStat(11);
C3->Print("gamma_Det2P.png");
*/



}
