/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 */

// This is needed for calling standalone classes
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <vector>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"

// FemtoDst headers
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoDstReader.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoDst.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoEvent.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoTrack.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

//_________________
void FemtoDstAnalyzerV2(const Char_t *inFile = "AuAu27GeV/AuAu27_ar.list", const Char_t *outFileName = "outFile27.root",const Char_t *flatFileName="Result27Recenter.root",const Char_t *resolutionFileName="Result27Recenter.root") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  //WEST
  TH1D *H_Psi2_west = new TH1D("H_Psi2_west","#psi_{2} of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_west = new TH1D("H_Psi3_west","#psi_{3} of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi2_recenter_west = new TH1D("H_Psi2_recenter_west","#psi^{recenter}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_recenter_west = new TH1D("H_Psi3_recenter_west","#psi^{recenter}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  TH1D *H_Psi2_flat_west = new TH1D("H_Psi2_flat_west","#psi^{rec+flat}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_flat_west = new TH1D("H_Psi3_flat_west","#psi^{rec+flat}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  //EAST
  TH1D *H_Psi2_east = new TH1D("H_Psi2_east","#psi_{2} of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_east = new TH1D("H_Psi3_east","#psi_{3} of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi2_recenter_east = new TH1D("H_Psi2_recenter_east","#psi^{recenter}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_recenter_east = new TH1D("H_Psi3_recenter_east","#psi^{recenter}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  TH1D *H_Psi2_flat_east = new TH1D("H_Psi2_flat_east","#psi^{rec+flat}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_flat_east = new TH1D("H_Psi3_flat_east","#psi^{rec+flat}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );

  //Profiles for v2 and v3
  TProfile2D *P_v2_pt = new TProfile2D("P_v2_pt","Profile of p_{t} versus v_{2})", 50, 0.2, 5.2, 9, -0.5, 8.5);
  TProfile2D *P_v3_pt = new TProfile2D("P_v3_pt","Profile of p_{t} versus v_{3}", 50, 0.2, 5.2, 9, -0.5, 8.5);
  TProfile *P_v2_cent = new TProfile("P_v2_cent","Profile of Centrality versus v_{2}", 9, -0.5, 8.5);
  TProfile *P_v3_cent = new TProfile("P_v3_cent","Profile of Centrality versus v_{3}", 9, -0.5, 8.5);

  //----------READING FLATENING-------------------------
  TFile *readf = new TFile(flatFileName,"READ");
  //_west
  TProfile2D *P3_Qx2_west = (TProfile2D*)readf->Get("P2_Qx2_cent_RunID_west");
  TProfile2D *P3_Qy2_west = (TProfile2D*)readf->Get("P2_Qy2_cent_RunID_west");
  TProfile2D *P3_Qx3_west = (TProfile2D*)readf->Get("P2_Qx3_cent_RunID_west");
  TProfile2D *P3_Qy3_west = (TProfile2D*)readf->Get("P2_Qy3_cent_RunID_west");
  //_east
  TProfile2D *P3_Qx2_east = (TProfile2D*)readf->Get("P2_Qx2_cent_RunID_east");
  TProfile2D *P3_Qy2_east = (TProfile2D*)readf->Get("P2_Qy2_cent_RunID_east");
  TProfile2D *P3_Qx3_east = (TProfile2D*)readf->Get("P2_Qx3_cent_RunID_east");
  TProfile2D *P3_Qy3_east = (TProfile2D*)readf->Get("P2_Qy3_cent_RunID_east");

  //WEST
  // Profile for cos and sin at v2
  TProfile2D *p_sin_v2_west[4];
  TProfile2D *p_cos_v2_west[4];
  for (int i=0; i<4; i++){
    p_sin_v2_west[i] = (TProfile2D*)readf -> Get(Form("sin_v2_west_prof_k=%i",i+1));
  }
  for (int i=0; i<4; i++){
    p_cos_v2_west[i] = (TProfile2D*)readf -> Get(Form("cos_v2_west_prof_k=%i",i+1));
  }
  // Profile for cos and sin at v3
  TProfile2D *p_sin_v3_west[4];
  TProfile2D *p_cos_v3_west[4];
  for (int i=0; i<4; i++){
    p_sin_v3_west[i] = (TProfile2D*)readf -> Get(Form("sin_v3_west_prof_k=%i",i+1));
  }
  for (int i=0; i<4; i++){
    p_cos_v3_west[i] = (TProfile2D*)readf -> Get(Form("cos_v3_west_prof_k=%i",i+1));
  }  

  //EAST
  // Profile for cos and sin at v2
  TProfile2D *p_sin_v2_east[4];
  TProfile2D *p_cos_v2_east[4];
  for (int i=0; i<4; i++){
    p_sin_v2_east[i] = (TProfile2D*)readf -> Get(Form("sin_v2_east_prof_k=%i",i+1));
  }
  for (int i=0; i<4; i++){
    p_cos_v2_east[i] = (TProfile2D*)readf -> Get(Form("cos_v2_east_prof_k=%i",i+1));
  }
  // Profile for cos and sin at v3
  TProfile2D *p_sin_v3_east[4];
  TProfile2D *p_cos_v3_east[4];
  for (int i=0; i<4; i++){
    p_sin_v3_east[i] = (TProfile2D*)readf -> Get(Form("sin_v3_east_prof_k=%i",i+1));
  }
  for (int i=0; i<4; i++){
    p_cos_v3_east[i] = (TProfile2D*)readf -> Get(Form("cos_v3_east_prof_k=%i",i+1));
  }  


  //----------READING RESOLUTION-------------------------  
  TFile *readf2 = new TFile(resolutionFileName,"READ");
  TProfile *P2_square_resolution2 = (TProfile*)readf2 -> Get("P_square_resolution2");
  TProfile *P2_square_resolution3 = (TProfile*)readf2 -> Get("P_square_resolution3");

  // //----------READING RESOLUTION-------------------------  
  // TFile *readf2 = new TFile(resolutionFileName,"READ");
  //TProfile *P2_square_resolution2 = (TProfile*)readf -> Get("P_square_resolution2");
  //TProfile *P2_square_resolution3 = (TProfile*)readf -> Get("P_square_resolution3");



  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  Double_t omega=0, Qx2=0, Qy2=0, Qx2_recenter=0, Qy2_recenter=0, Psi2=0, Psi2_recenter=0, delta_Psi2=0, Psi2_flat=0;
  Double_t Qx3=0, Qy3=0, Qx3_recenter=0, Qy3_recenter=0, Psi3=0, Psi3_recenter=0, delta_Psi3=0, Psi3_flat=0;
  Double_t Q2_west=0, Qx2_west=0, Qy2_west=0, Psi2_west=0, Q3_west=0, Qx3_west=0, Qy3_west=0, Psi3_west=0, Qx2_recenter_west=0, Qy2_recenter_west=0, Qx3_recenter_west=0, Qy3_recenter_west=0, Psi2_recenter_west=0, Psi3_recenter_west=0, Psi2_flat_west=0, Psi3_flat_west=0, delta_Psi2_west=0, delta_Psi3_west=0;
  Double_t Q2_east=0, Qx2_east=0, Qy2_east=0, Psi2_east=0, Q3_east=0, Qx3_east=0, Qy3_east=0, Psi3_east=0, Qx2_recenter_east=0, Qy2_recenter_east=0, Qx3_recenter_east=0, Qy3_recenter_east=0, Psi3_recenter_east=0, Psi2_recenter_east=0, Psi2_flat_east=0, Psi3_flat_east=0, delta_Psi2_east=0, delta_Psi3_east=0;
  Double_t N_west=0, N_east=0;
  Double_t R2_cos = 0, R3_cos = 0, v2 = 0, v3 = 0, weight = 0;;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    //std::cout << "Working on event #[" << (iEvent+1)
    //  	      << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
    StFemtoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    // Return primary vertex position
    TVector3 pVtx = event->primaryVertex();

    // Reject vertices that are far from the central membrane along the beam
    if( TMath::Abs( pVtx.Z() ) > 70. ) continue;
    if( TMath::Abs( pow(pVtx.X(), 2)+ pow(pVtx.Y(), 2)) > 2. ) continue;
    //if (event->vpdVz() == 0.0) continue;

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();

    N_west = 0;
    N_east = 0;
    v2 = 0;
    v3 = 0;

    //Q-vector cleaning in new event
    Qx2=0; Qy2=0; Qx2_recenter=0; Qy2_recenter=0; Psi2=0; Psi2_recenter=0; Psi2_flat=0; delta_Psi2=0;
    Qx3=0; Qy3=0; Qx3_recenter=0; Qy3_recenter=0; Psi3=0; Psi3_recenter=0; Psi3_flat=0; delta_Psi3=0;

    Qx2_west=0; Qy2_west=0; Qx2_recenter_west=0; Qy2_recenter_west=0; Psi2_west=0; Psi2_recenter_west=0; Psi2_flat_west=0; delta_Psi2_west=0;
    Qx3_east=0; Qy3_east=0; Qx3_recenter_east=0; Qy3_recenter_east=0; Psi3_east=0; Psi3_recenter_east=0; Psi3_flat_east=0; delta_Psi3_east=0;

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;

      // Must be a primary track
      if ( !femtoTrack->isPrimary() ) continue;
      if( (femtoTrack->dEdx()) == 0 ) continue;
      // Simple single-track cut
      if( femtoTrack->gMom().Mag() < 0.1 || femtoTrack->gDCA(pVtx).Mag() > 2. ) {
        continue;
      }
      if( femtoTrack -> p() < 0.15 || femtoTrack -> p() > 5 || 
      TMath::Abs( femtoTrack -> eta() ) > 1 || femtoTrack -> nHits() < 15 ) { /**/ //15 из 45 падов сработали, при eta>1 эффективность сильно падает
        continue;
      }

      if(femtoTrack->pt()>2.0 || femtoTrack->pt()<0.2){ 
        continue;
      }
      
      omega=femtoTrack->pt();
      if(femtoTrack->eta()>0.05){
        Qx2_east+= omega*cos(2*(femtoTrack->phi()) );
        Qy2_east+= omega*sin(2*(femtoTrack->phi()) );
        Qx3_east+= omega*cos(3*(femtoTrack->phi()) );
        Qy3_east+= omega*sin(3*(femtoTrack->phi()) );
        N_east++;
      }

      if(femtoTrack->eta()<-0.05){
        Qx2_west+= omega*cos(2*(femtoTrack->phi()) );
        Qy2_west+= omega*sin(2*(femtoTrack->phi()) );
        Qx3_west+= omega*cos(3*(femtoTrack->phi()) );
        Qy3_west+= omega*sin(3*(femtoTrack->phi()) );
        N_west++;
      }


      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {

      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    if( N_west != 0 ){
      Qx2_west = Qx2_west / N_west; 
      Qy2_west = Qy2_west / N_west;
      Qx3_west = Qx3_west / N_west;
      Qy3_west = Qy3_west / N_west;
    }
    if( N_east != 0 ){
      Qx2_east = Qx2_east / N_east; 
      Qy2_east = Qy2_east / N_east;
      Qx3_east = Qx3_east / N_east;
      Qy3_east = Qy3_east / N_east;
    }

  //Finding Q2,Q3,Psi2,Psi3
  //------west arm Q
  if( N_west != 0 ){
    Psi2_west=1.0/2.0*(TMath::ATan2(Qy2_west, Qx2_west)); //(TMath::Sqrt(event->refMult()))
    Psi3_west=1.0/3.0*(TMath::ATan2(Qy3_west, Qx3_west)); 
    Q2_west=TMath::Sqrt(pow(Qx2_west,2)+pow(Qy2_west,2));
    Q3_west=TMath::Sqrt(pow(Qx3_west,2)+pow(Qy3_west,2));
  }
  //-----east Q
  if( N_east != 0 ){
    Psi2_east=1.0/2.0*(TMath::ATan2(Qy2_east, Qx2_east)); //(TMath::Sqrt(event->refMult()))
    Psi3_east=1.0/3.0*(TMath::ATan2(Qy3_east, Qx3_east)); 
    Q2_east=TMath::Sqrt(pow(Qx2_east,2)+pow(Qy2_east,2));
    Q3_east=TMath::Sqrt(pow(Qx3_east,2)+pow(Qy3_east,2));
  }
  //----_west
  if( N_west != 0 ){
    Qx2_recenter_west=Qx2_west- P3_Qx2_west->GetBinContent(P3_Qx2_west->FindBin(event->runId(), event->cent9()));
    Qy2_recenter_west=Qy2_west- P3_Qy2_west->GetBinContent(P3_Qy2_west->FindBin(event->runId(), event->cent9()));
    Qx3_recenter_west=Qx3_west- P3_Qx3_west->GetBinContent(P3_Qx3_west->FindBin(event->runId(), event->cent9()));
    Qy3_recenter_west=Qy3_west- P3_Qy3_west->GetBinContent(P3_Qy3_west->FindBin(event->runId(), event->cent9()));
  }
  //----east
  if( N_east != 0 ){
    Qx2_recenter_east=Qx2_east- P3_Qx2_east->GetBinContent(P3_Qx2_east->FindBin(event->runId(), event->cent9()));
    Qy2_recenter_east=Qy2_east- P3_Qy2_east->GetBinContent(P3_Qy2_east->FindBin(event->runId(), event->cent9()));
    Qx3_recenter_east=Qx3_east- P3_Qx3_east->GetBinContent(P3_Qx3_east->FindBin(event->runId(), event->cent9()));
    Qy3_recenter_east=Qy3_east- P3_Qy3_east->GetBinContent(P3_Qy3_east->FindBin(event->runId(), event->cent9()));
  }
  //_west
  if( N_west != 0 ){
    Psi2_recenter_west=1.0/2.0*(TMath::ATan2(Qy2_recenter_west, Qx2_recenter_west)); //(TMath::Sqrt(event->refMult()))
    Psi3_recenter_west=1.0/3.0*(TMath::ATan2(Qy3_recenter_west, Qx3_recenter_west));
  }
  //_east
  if( N_east != 0 ){
    Psi2_recenter_east=1.0/2.0*(TMath::ATan2(Qy2_recenter_east, Qx2_recenter_east)); //(TMath::Sqrt(event->refMult()))
    Psi3_recenter_east=1.0/3.0*(TMath::ATan2(Qy3_recenter_east, Qx3_recenter_east));
  }
  //WEST
  if( N_west != 0 ){
    for (int i=0; i<4; i++){
        delta_Psi2_west += 2.0/( (Double_t)i+1.0)*(p_cos_v2_west[i]->GetBinContent(p_cos_v2_west[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*2*Psi2_recenter_west)
        -p_sin_v2_west[i]->GetBinContent(p_sin_v2_west[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*2*Psi2_recenter_west) );
      }

      Psi2_flat_west = Psi2_recenter_west + 1.0/2.0*delta_Psi2_west;

    for (int i=0; i<4; i++){
      delta_Psi3_west += 2.0/( (Double_t)i+1.0)*(p_cos_v3_west[i]->GetBinContent(p_cos_v3_west[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*3*Psi3_recenter_west)
      -p_sin_v3_west[i]->GetBinContent(p_sin_v3_west[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*3*Psi3_recenter_west) );
    }

      Psi3_flat_west = Psi3_recenter_west +1.0/3.0*delta_Psi3_west;
  }

  //EAST
  if( N_east != 0 ){
    for (int i=0; i<4; i++){
        delta_Psi2_east += 2.0/( (Double_t)i+1.0)*(p_cos_v2_east[i]->GetBinContent(p_cos_v2_east[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*2*Psi2_recenter_east)
        -p_sin_v2_east[i]->GetBinContent(p_sin_v2_east[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*2*Psi2_recenter_east) );
      }

      Psi2_flat_east = Psi2_recenter_east + 1.0/2.0*delta_Psi2_east;

    for (int i=0; i<4; i++){
      delta_Psi3_east += 2.0/( (Double_t)i+1.0)*(p_cos_v3_east[i]->GetBinContent(p_cos_v3_east[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*3*Psi3_recenter_east)
      -p_sin_v3_east[i]->GetBinContent(p_sin_v3_east[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*3*Psi3_recenter_east) );
    }

      Psi3_flat_east = Psi3_recenter_east + 1.0/3.0*delta_Psi3_east;
  }


  // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if (!femtoTrack) continue;
      if ( !femtoTrack->isPrimary() ) continue;
      if( (femtoTrack->dEdx()) == 0 ) continue;
      if( femtoTrack->gMom().Mag() < 0.1 || femtoTrack->gDCA(pVtx).Mag() > 2. ) {
        continue;
      }
      if( femtoTrack -> p() < 0.15 || femtoTrack -> p() > 5 || TMath::Abs( femtoTrack -> eta() ) > 1 || femtoTrack -> nHits() < 15 ) {
        continue;
      }
      if (femtoTrack->pt()>0 && femtoTrack->pt()<2){
        omega = femtoTrack->pt();}
      else {
        omega = 2;
      }
      if(femtoTrack->eta()>0.05){
        v2 = TMath::Cos(2*(femtoTrack->phi()-Psi2_flat_east));
        v3 = TMath::Cos(3*(femtoTrack->phi()-Psi3_flat_east));
        N_east++;
      }

      if(femtoTrack->eta()<-0.05){
        v2 = TMath::Cos(2*(femtoTrack->phi()-Psi2_flat_west));
        v3 = TMath::Cos(3*(femtoTrack->phi()-Psi3_flat_west));
        N_west++;
      }
      P_v2_pt -> Fill(femtoTrack->pt(), event->cent9(), v2/TMath::Sqrt(P2_square_resolution2 -> GetBinContent(P2_square_resolution2->FindBin(event->cent9()))), omega);
      P_v3_pt -> Fill(femtoTrack->pt(), event->cent9(), v3/TMath::Sqrt(P2_square_resolution3 -> GetBinContent(P2_square_resolution3->FindBin(event->cent9()))), omega);
      P_v2_cent -> Fill(event->cent9(), v2/TMath::Sqrt(P2_square_resolution2 -> GetBinContent(P2_square_resolution2->FindBin(event->cent9()))), omega);
      P_v3_cent -> Fill(event->cent9(), v3/TMath::Sqrt(P2_square_resolution3 -> GetBinContent(P2_square_resolution3->FindBin(event->cent9()))), omega);
    }

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  //Saving our files
  TFile *savefile = new TFile(outFileName, "RECREATE");

  //_west
  H_Psi2_west->Write();
  H_Psi3_west->Write();
  H_Psi2_recenter_west->Write();
  H_Psi3_recenter_west->Write();
  H_Psi2_flat_west->Write();
  H_Psi3_flat_west->Write();
  //_east
  H_Psi2_east->Write();
  H_Psi3_east->Write();
  H_Psi2_recenter_east->Write();
  H_Psi3_recenter_east->Write();
  H_Psi2_flat_east->Write();
  H_Psi3_flat_east->Write();

  P_v2_pt->Write();
  P_v3_pt->Write();
  P_v2_cent->Write();
  P_v3_cent->Write();

  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}