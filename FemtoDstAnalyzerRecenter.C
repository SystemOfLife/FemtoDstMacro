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
void FemtoDstAnalyzerRecenter(const Char_t *inFile = "AuAu27GeV/AuAu27_ar.list", const Char_t *outFileName = "outFile27.root", const Char_t *recFileName="Result27Profile2D.root") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  TH1D *H_Qx2_recenter = new TH1D("H_Qx2_recenter", "Q^{recenter}_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy2_recenter = new TH1D("H_Qy2_recenter", "Q^{recenter}_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qx3_recenter = new TH1D("H_Qx3_recenter", "Q^{recenter}_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy3_recenter = new TH1D("H_Qy3_recenter", "Q^{recenter}_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Psi2 = new TH1D("H_Psi2","#psi_{2} of Event plane of 2 harmonic; #psi,[rad];N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3 = new TH1D("H_Psi3","#psi_{3} of Event plane of 3 harmonic; #psi,[rad];N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi2_recenter = new TH1D("H_Psi2_recenter","#psi^{recenter}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_recenter = new TH1D("H_Psi3_recenter","#psi^{recenter}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  TProfile2D *P2_Qx2_cent_RunID = new TProfile2D("P2_Qx2_cent_RunID","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P2_Qy2_cent_RunID = new TProfile2D("P2_Qy2_cent_RunID","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P2_Qx3_cent_RunID = new TProfile2D("P2_Qx3_cent_RunID","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P2_Qy3_cent_RunID = new TProfile2D("P2_Qy3_cent_RunID","Profile of Q_{3}_{y} on Centrality vs RunID", 0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );

  // Profiles for cos and sin at v2
  TProfile2D *sin_v2[4];
  for (int i=0; i<4; i++){
    sin_v2[i] = new TProfile2D(Form("sin_v2_prof_k=%i",i+1), Form("sin_v2_prof_k=%i",i+1), 0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }
  TProfile2D *cos_v2[4];
  for (int i=0; i<4; i++){
    cos_v2[i] = new TProfile2D(Form("cos_v2_prof_k=%i",i+1), Form("cos_v2_prof_k=%i",i+1), 0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }
  // Profiles for cos and sin at v3
  TProfile2D *sin_v3[4];
  for (int i=0; i<4; i++){
    sin_v3[i] = new TProfile2D(Form("sin_v3_prof_k=%i",i+1), Form("sin_v3_prof_k=%i",i+1), 0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }
  TProfile2D *cos_v3[4];
  for (int i=0; i<4; i++){
    cos_v3[i] = new TProfile2D(Form("cos_v3_prof_k=%i",i+1), Form("cos_v3_prof_k=%i",i+1), 0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }

  TFile *readf = new TFile(recFileName,"READ");
  TProfile2D *P_Qx2 = (TProfile2D*)readf->Get("P_Qx2_cent_RunID");
  TProfile2D *P_Qy2 = (TProfile2D*)readf->Get("P_Qy2_cent_RunID");
  TProfile2D *P_Qx3 = (TProfile2D*)readf->Get("P_Qx3_cent_RunID");
  TProfile2D *P_Qy3 = (TProfile2D*)readf->Get("P_Qy3_cent_RunID");

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

  Double_t omega=0, Qx2=0, Qy2=0, Qx2_recenter=0, Qy2_recenter=0, Psi2=0, Psi2_recenter=0;
  Double_t Qx3=0, Qy3=0, Qx3_recenter=0, Qy3_recenter=0, Psi3=0, Psi3_recenter=0;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    // std::cout << "Working on event #[" << (iEvent+1)
    //   	      << "/" << events2read << "]" << std::endl;

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
    if( TMath::Abs( pVtx.Z() ) > 40. ) continue;
    if( TMath::Abs( pow(pVtx.X(), 2)+ pow(pVtx.Y(), 2)) > 2. ) continue;
    if (event->vpdVz() == 0.0) continue;

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();

    //Q-vector cleaning in new event
    Qx2=0, Qy2=0, Qx2_recenter=0, Qy2_recenter=0, Psi2=0, Psi2_recenter=0;
    Qx3=0, Qy3=0, Qx3_recenter=0, Qy3_recenter=0, Psi3=0, Psi3_recenter=0;

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;

      // Must be a primary track
      if ( !femtoTrack->isPrimary() ) continue;

      if( (femtoTrack->dEdx()) == 0 ) continue;

      // Simple single-track cut
      if( femtoTrack->gMom().Mag() < 0.1 || femtoTrack->gDCA(pVtx).Mag() > 3. ) {
        continue;
      }
      if( femtoTrack -> p() < 0.1 || femtoTrack -> p() > 10 || TMath::Abs( femtoTrack -> eta() ) > 1 || femtoTrack -> nHits() < 15 ) { /**/ //15 из 45 падов сработали, при eta>1 эффективность сильно падает
        continue;
      }

      if(femtoTrack->pt()>2.0){ 
        omega=2.0;
      }
      else { 
        omega=femtoTrack->pt();
      }
      //Qx and Qy
      Qx2+= omega*cos(2*(femtoTrack->phi()) );
      Qy2+= omega*sin(2*(femtoTrack->phi()) );
      Qx3+= omega*cos(3*(femtoTrack->phi()) );
      Qy3+= omega*sin(3*(femtoTrack->phi()) );

      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {

      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    //Finding Q2,Q3,Psi2,Psi3
  Psi2=1.0/2.0*(TMath::ATan2(Qy2, Qx2)); //(TMath::Sqrt(event->refMult()))
  Psi3=1.0/3.0*(TMath::ATan2(Qy3, Qx3)); 

  Qx2_recenter=Qx2- P_Qx2->GetBinContent(P_Qx2->FindBin(event->runId(), event->cent9()));
  Qy2_recenter=Qy2- P_Qy2->GetBinContent(P_Qy2->FindBin(event->runId(), event->cent9()));
  Qx3_recenter=Qx3- P_Qx3->GetBinContent(P_Qx3->FindBin(event->runId(), event->cent9()));
  Qy3_recenter=Qy3- P_Qy3->GetBinContent(P_Qy3->FindBin(event->runId(), event->cent9()));

  Psi2_recenter=1.0/2.0*(TMath::ATan2(Qy2_recenter, Qx2_recenter)); //(TMath::Sqrt(event->refMult()))
  Psi3_recenter=1.0/3.0*(TMath::ATan2(Qy3_recenter, Qx3_recenter));

  //Filling profiles of cos and sin
  for (int i=0; i<4; i++){
    sin_v2[i] -> Fill(event->runId(), event->cent9(), TMath::Sin( (i+1)*2*Psi2_recenter ), 1);
  }    
  for (int i=0; i<4; i++){
    cos_v2[i] -> Fill(event->runId(), event->cent9(), TMath::Cos( (i+1)*2*Psi2_recenter ), 1);
  }
  for (int i=0; i<4; i++){
    sin_v3[i] -> Fill(event->runId(), event->cent9(), TMath::Sin( (i+1)*3*Psi3_recenter ), 1);
  }   
  for (int i=0; i<4; i++){
    cos_v3[i] -> Fill(event->runId(), event->cent9(), TMath::Cos( (i+1)*3*Psi3_recenter ), 1);
  }

  P2_Qx2_cent_RunID->Fill(event->runId(), event->cent9(), Qx2, 1);
  P2_Qy2_cent_RunID->Fill(event->runId(), event->cent9(), Qy2, 1);
  P2_Qx3_cent_RunID->Fill(event->runId(), event->cent9(), Qx3, 1);
  P2_Qy3_cent_RunID->Fill(event->runId(), event->cent9(), Qy3, 1);

  

  H_Qx2_recenter->Fill(Qx2_recenter);
  H_Qy2_recenter->Fill(Qy2_recenter);
  H_Qx3_recenter->Fill(Qx3_recenter);
  H_Qy3_recenter->Fill(Qy3_recenter);
  H_Psi2->Fill(Psi2);
  H_Psi3->Fill(Psi3);
  H_Psi2_recenter->Fill(Psi2_recenter);
  H_Psi3_recenter->Fill(Psi3_recenter);

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  //Saving our files
  TFile *savefile = new TFile(outFileName, "RECREATE");

  H_Qx2_recenter->Write();
  H_Qy2_recenter->Write();
  H_Qx3_recenter->Write();
  H_Qy3_recenter->Write();
  H_Psi2->Write();
  H_Psi3->Write();
  H_Psi2_recenter->Write();
  H_Psi3_recenter->Write();
  P2_Qx2_cent_RunID->Write();
  P2_Qy2_cent_RunID->Write();
  P2_Qx3_cent_RunID->Write();
  P2_Qy3_cent_RunID->Write();
  for (int i=0; i<4; i++){
    sin_v2[i]->Write();
    cos_v2[i]->Write();
    sin_v3[i]->Write();
    cos_v3[i]->Write();
  }

  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}
