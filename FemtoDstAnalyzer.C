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
R__LOAD_LIBRARY(/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/libStFemtoDst)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

//_________________
void FemtoDstAnalyzer(const Char_t *inFile = "AuAu27GeV/AuAu27_ar.list", const Char_t *outFileName = "outFile27.root") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/libStFemtoDst.so");
  #endif

  //TH1D, basic distribution:
  TH1D *H_massSqr = new TH1D ("H_massSqr","Squared mass; m^{2},[(GeV/c)^{2}]", 500, 0, 1.2 );
  TH2D *H_VxVy = new TH2D("H_VxVy","X vs Y vertex position; x[cm]; y[cm]", 1000, -2.5, 2.5, 1000, -2.5, 2.5);
  TH1D *H_invBeta = new TH1D ("H_invBeta", "1/beta", 500, -0.2, 4);
  TH1D *H_refMult = new TH1D ("H_refMult", "Multiplicity; N_{ch}; N_{count}",500, 0, 450);
  TH1D *H_cent9 = new TH1D ("H_cent9", "Centrality; Centrality bin; N_{count}", 40, 0, 20);
  TH1D *H_nHits = new TH1D ("H_nHits", "N_{hits} distribution ; N_{hits}; N_{count}", 100, 0, 45);
  TH1D *H_VtxZ = new TH1D("H_VtxZ","Vertex Z distribution; VtxZ,[cm]; N_{count}" ,500, -55, 55);
  TH1D *H_dEdx = new TH1D("H_dEdx","Energy loss; dE/dx,[a.u.]",500, 0, 0.00007);
  TH1D *H_Pt = new TH1D("H_Pt","Transverse momentum distribution; Pt,[GeV/c]", 1000 , 0, 6);
  TH1D *H_Px = new TH1D("H_Px","Px distribution; Px,[GeV/c]", 1000, 0, 6);
  TH1D *H_Py = new TH1D("H_Py","Py distribution; Py,[GeV/c]", 1000, 0,6);
  TH1D *H_Phi = new TH1D("H_Phi","Phi distribution; #phi [rad]", 1000, -3.5, 3.5);
  TH1D *H_Eta = new TH1D("H_Eta","Eta distribution; #eta, [a.u.]", 1000, -2, 2);
  TH1D *H_runID = new TH1D("H_runID","RunID distribution; Number of Run; N_{count}", 500, 10e6, 13e6);
  //Q-vectors hystograms:
  TH1D *H_Qx2 = new TH1D("H_Qx2","Q_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy2 = new TH1D("H_Qy2","Q_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q2 = new TH1D("H_Q2","Q_{2}-vector; Q_{2},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qx3 = new TH1D("H_Qx3","Q_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy3 = new TH1D("H_Qy3","Q_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q3 = new TH1D("H_Q3","Q_{3}-vector; Q_{3},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Psi2 = new TH1D("H_Psi2","#psi of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3 = new TH1D("H_Psi3","#psi of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
   //Q-vectors WEST hystograms:
  TH1D *H_Qx2_west = new TH1D("H_Qx2_west","Q_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy2_west = new TH1D("H_Qy2_west","Q_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q2_west = new TH1D("H_Q2_west","Q_{2}-vector; Q_{2},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qx3_west = new TH1D("H_Qx3_west","Q_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy3_west = new TH1D("H_Qy3_west","Q_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q3_west = new TH1D("H_Q3_west","Q_{3}-vector; Q_{3},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Psi2_west = new TH1D("H_Psi2_west","#psi of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_west = new TH1D("H_Psi3_west","#psi of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  //Q-vectors EAST hystograms:
  TH1D *H_Qx2_east = new TH1D("H_Qx2_east","Q_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy2_east = new TH1D("H_Qy2_east","Q_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q2_east = new TH1D("H_Q2_east","Q_{2}-vector; Q_{2},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qx3_east = new TH1D("H_Qx3_east","Q_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Qy3_east = new TH1D("H_Qy3_east","Q_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q3_east = new TH1D("H_Q3_east","Q_{3}-vector; Q_{3},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Psi2_east = new TH1D("H_Psi2_east","#psi of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_east = new TH1D("H_Psi3_east","#psi of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );

  //Invariant mass distribution:
  // TH1D *H_MinvPiTPC = new TH1D ("H_MinvPiTPC", "Invariant mass of #pi^{+-}; M_{inv},[GeV/c^{2}]; N_{counts}", 500, 0.25, 1.2);
  // TH1D *H_MinvKaTPC = new TH1D ("H_MinvKaTPC", "Invariant mass of K^{+}K^{-}; M_{inv},[GeV/c^{2}]; N_{counts}", 500, 0.98, 1.4);
  // TH1D *H_MinvPPiPosTPC = new TH1D ("H_MinvPPiPosTPC", "Invariant mass of Proton,#pi^{+}; M_{inv},[GeV/c^{2}]; N_{counts}", 500, 0.9, 1.6);
  // TH1D *H_MinvPiTOF = new TH1D ("H_MinvPiTOF", "Invariant mass of #pi^{+-}; M_{inv},[GeV/c^{2}]; N_{counts}", 500, 0.25, 1.2);
  // TH1D *H_MinvKaTOF = new TH1D ("H_MinvKaTOF", "Invariant mass of K^{+}K^{-}; M_{inv} [GeV/c^{2}]; N_{counts}", 500, 0.98, 1.4);
  // TH1D *H_MinvPPiPosTOF = new TH1D ("H_MinvPPiPosTOF", "Invariant mass of Proton,#pi^{+}; M_{inv},[GeV/c^{2}]; N_{counts}", 500, 0.9, 1.6);

  //TH2F
  TH2D *H_PtEta = new TH2D("H_PtEta","P_{t} vs Eta; P_{t},[GeV/c]; #eta", 1000, 0,6,   1000,-2,2);
  TH2D *H_PtPhi = new TH2D("H_PtPhi","P_{t} vs Phi; P_{t},[GeV/c] ; #phi[rad]", 1000,0,6,  1000,-3.2,3.2);
  TH2D *H_EtaPhi = new TH2D("H_EtaPhi","Eta vs Phi; #eta; #phi[rad]", 1000,-3.2,3.2,  1000,-2,2);
  TH2D *H_dEdxQP = new TH2D ("H_dEdxQP","dE/dx vs Q*P; Q*P,[GeV/c]; dE/dx,[a.u.]", 1000,-3.5,3.5,   1000,0,0.00003);
  TH2D *H_m2QP = new TH2D ("H_m2QP","m^{2} vs Q*P; Q*P,[GeV/c]; m^{2},[(GeV/c)^{2}]", 500,-2.5,2.5,  500,-1,1.3);
  TH2D *H_nSigK_QP = new TH2D ("H_nSigK_QP", "n#sigma Kaon vs Q*P; Q*P,[GeV/c]; n#sigma Kaon,[a.u.]", 500,-2.5,2.5,   500,-8,8);
  TH2D *H_nSigPi_QP = new TH2D ("H_nSigPi_QP", "n#sigma Pion vs Q*P; Q*P,[GeV/c]; n#sigma Pion,[a.u.]", 500,-2.5,2.5,   500,-8,8);
  TH2D *H_nSigP_QP = new TH2D ("H_nSigP_QP", "n#sigma Proton vs Q*P; Q*P,[GeV/c]; n#sigma Proton,[a.u.]", 500,-2.5,2.5,   500,-8,8);

  //TProfiles
  TProfile *P_Q2_cent = new TProfile("P_Q2_cent","Centrality vs Q_{2}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx2_cent = new TProfile("P_Qx2_cent","Centrality vs Q_{2}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy2_cent = new TProfile("P_Qy2_cent","Centrality vs Q_{2}_{y}",250,-0.5,10,  -100,100 );
  TProfile *P_Q3_cent = new TProfile("P_Q3_cent","Centrality vs Q_{3}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx3_cent = new TProfile("P_Qx3_cent","Centrality vs Q_{3}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy3_cent = new TProfile("P_Qy3_cent","Centrality vs Q_{3}_{y}",250,-0.5,10,  -100,100 );
  TProfile2D *P_Qx2_cent_RunID = new TProfile2D("P_Qx2_cent_RunID","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy2_cent_RunID = new TProfile2D("P_Qy2_cent_RunID","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qx3_cent_RunID = new TProfile2D("P_Qx3_cent_RunID","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy3_cent_RunID = new TProfile2D("P_Qy3_cent_RunID","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  //-----west
  TProfile *P_Q2_cent_west = new TProfile("P_Q2_cent_west","Centrality vs Q_{2}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx2_cent_west = new TProfile("P_Qx2_cent_west","Centrality vs Q_{2}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy2_cent_west = new TProfile("P_Qy2_cent_west","Centrality vs Q_{2}_{y}",250,-0.5,10,  -100,100 );
  TProfile *P_Q3_cent_west = new TProfile("P_Q3_cent_west","Centrality vs Q_{3}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx3_cent_west = new TProfile("P_Qx3_cent_west","Centrality vs Q_{3}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy3_cent_west = new TProfile("P_Qy3_cent_west","Centrality vs Q_{3}_{y}",250,-0.5,10,  -100,100 );
  TProfile2D *P_Qx2_cent_RunID_west = new TProfile2D("P_Qx2_cent_RunID_west","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy2_cent_RunID_west = new TProfile2D("P_Qy2_cent_RunID_west","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qx3_cent_RunID_west = new TProfile2D("P_Qx3_cent_RunID_west","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy3_cent_RunID_west = new TProfile2D("P_Qy3_cent_RunID_west","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  //-----east
  TProfile *P_Q2_cent_east= new TProfile("P_Q2_cent_east","Centrality vs Q_{2}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx2_cent_east = new TProfile("P_Qx2_cent_east","Centrality vs Q_{2}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy2_cent_east = new TProfile("P_Qy2_cent_east","Centrality vs Q_{2}_{y}",250,-0.5,10,  -100,100 );
  TProfile *P_Q3_cent_east= new TProfile("P_Q3_cent_east","Centrality vs Q_{3}",50,-0.5,10,  -100,100 );
  TProfile *P_Qx3_cent_east = new TProfile("P_Qx3_cent_east","Centrality vs Q_{3}_{x}",250,-0.5,10,  -100,100 );
  TProfile *P_Qy3_cent_east = new TProfile("P_Qy3_cent_east","Centrality vs Q_{3}_{y}",250,-0.5,10,  -100,100 );
  TProfile2D *P_Qx2_cent_RunID_east = new TProfile2D("P_Qx2_cent_RunID_east","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy2_cent_RunID_east = new TProfile2D("P_Qy2_cent_RunID_east","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qx3_cent_RunID_east = new TProfile2D("P_Qx3_cent_RunID_east","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy3_cent_RunID_east = new TProfile2D("P_Qy3_cent_RunID_east","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  

  Double_t omega=0,Q2=0,Qx2=0,Qy2=0,Psi2=0;
  Double_t Q3=0,Qx3=0,Qy3=0,Psi3=0;
  Double_t Q2_west=0, Qx2_west=0, Qy2_west=0, Psi2_west=0, Q3_west=0, Qx3_west=0, Qy3_west=0, Psi3_west=0;
  Double_t Q2_east=0, Qx2_east=0, Qy2_east=0, Psi2_east=0, Q3_east=0, Qx3_east=0, Qy3_east=0, Psi3_east=0;
  Double_t N_east=0, N_west=0;


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

  

  // Loop over events/////////////////////////////////////////////////
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    //comment section below
    //
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
    if( TMath::Abs( pVtx.Z() ) > 70. ) continue;
    if( TMath::Abs( pow(pVtx.X(), 2)+ pow(pVtx.Y(), 2)) > 2. ) continue;
    //if (event->vpdVz() == 0.0) continue;

    /////////
    //Filling event Histograms
    H_cent9->Fill(event->cent9());
    H_VtxZ->Fill(event-> vpdVz());
    H_refMult->Fill(event->refMult());
    H_VxVy->Fill(pVtx(0),pVtx(1));
    H_runID->Fill(event->runId());
    ////////////

    //Q-vector cleaning in new event
    Q2=0; Qx2=0; Qy2=0; Psi2=0; Q3=0; Qx3=0; Qy3=0; Psi3=0;
    Q2_west=0; Qx2_west=0; Qy2_west=0; Psi2_west=0; Q3_west=0; Qx3_west=0; Qy3_west=0; Psi3_west=0;
    Q2_east=0; Qx2_east=0; Qy2_east=0; Psi2_east=0; Q3_east=0; Qx3_east=0; Qy3_east=0; Psi3_east=0;
    N_east= 0; N_west=0;



    // Track analysis
    Int_t nTracks = dst->numberOfTracks();



    // Track loop///////////////////////////////
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
      //Qx and Qy
      Qx2+= omega*cos(2*(femtoTrack->phi()) );
      Qy2+= omega*sin(2*(femtoTrack->phi()) );
      Qx3+= omega*cos(3*(femtoTrack->phi()) );
      Qy3+= omega*sin(3*(femtoTrack->phi()) );
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

      //Track histograms
      H_dEdx->Fill(femtoTrack->dEdx());
      H_Pt->Fill(femtoTrack->pt());
      H_Px->Fill((femtoTrack->pt())*cos(femtoTrack->phi()));
      H_Py->Fill((femtoTrack->pt())*sin(femtoTrack->phi()));
      H_Phi->Fill(femtoTrack->phi());
      H_Eta->Fill(femtoTrack->eta());
      H_PtEta->Fill(femtoTrack->pt(), femtoTrack->eta());
      H_PtPhi->Fill(femtoTrack->pt(), femtoTrack->phi());
      H_EtaPhi->Fill(femtoTrack->phi(), femtoTrack->eta());
      H_nHits->Fill(femtoTrack->nHits()); // начальные распределения
      H_dEdxQP->Fill(femtoTrack->pt()*femtoTrack->charge(), femtoTrack->dEdx());
      H_nSigK_QP->Fill(femtoTrack->pt()*femtoTrack->charge(), femtoTrack->nSigmaKaon());
      H_nSigPi_QP->Fill(femtoTrack->pt()*femtoTrack->charge(), femtoTrack->nSigmaPion());
      H_nSigP_QP->Fill(femtoTrack->pt()*femtoTrack->charge(), femtoTrack->nSigmaProton());

      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {
        H_massSqr->Fill(femtoTrack->massSqr());
        H_m2QP->Fill(femtoTrack->p()*femtoTrack->charge(), femtoTrack->massSqr());
        H_invBeta->Fill(femtoTrack->invBeta());

      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)


    //Finding Q2,Q3,Psi2,Psi3
    Psi2=1.0/2.0*(TMath::ATan2(Qy2, Qx2)); //(TMath::Sqrt(event->refMult()))
    Psi3=1.0/3.0*(TMath::ATan2(Qy3, Qx3)); 
    Q2=TMath::Sqrt(pow(Qx2,2)+pow(Qy2,2));
    Q3=TMath::Sqrt(pow(Qx3,2)+pow(Qy3,2));
    //------west arm Q
    Psi2_west=1.0/2.0*(TMath::ATan2(Qy2_west, Qx2_west)); //(TMath::Sqrt(event->refMult()))
    Psi3_west=1.0/3.0*(TMath::ATan2(Qy3_west, Qx3_west)); 
    Q2_west=TMath::Sqrt(pow(Qx2_west,2)+pow(Qy2_west,2));
    Q3_west=TMath::Sqrt(pow(Qx3_west,2)+pow(Qy3_west,2));
    //-----east Q
    Psi2_east=1.0/2.0*(TMath::ATan2(Qy2_east, Qx2_east)); //(TMath::Sqrt(event->refMult()))
    Psi3_east=1.0/3.0*(TMath::ATan2(Qy3_east, Qx3_east)); 
    Q2_east=TMath::Sqrt(pow(Qx2_east,2)+pow(Qy2_east,2));
    Q3_east=TMath::Sqrt(pow(Qx3_east,2)+pow(Qy3_east,2));

    //Fillig Profiles
    P_Q2_cent->Fill(event->cent9(), Q2);
    P_Qx2_cent->Fill(event->cent9(), Qx2);
    P_Qy2_cent->Fill(event->cent9(), Qy2);
    P_Q3_cent->Fill(event->cent9(), Q3);
    P_Qx3_cent->Fill(event->cent9(), Qx3);
    P_Qy3_cent->Fill(event->cent9(), Qy3);
    P_Qx2_cent_RunID->Fill(event->runId(), event->cent9(), Qx2, 1);
    P_Qy2_cent_RunID->Fill(event->runId(), event->cent9(), Qy2, 1);
    P_Qx3_cent_RunID->Fill(event->runId(), event->cent9(), Qx3, 1);
    P_Qy3_cent_RunID->Fill(event->runId(), event->cent9(), Qy3, 1);
    //----west
    P_Q2_cent_west->Fill(event->cent9(), Q2_west);
    P_Qx2_cent_west->Fill(event->cent9(), Qx2_west);
    P_Qy2_cent_west->Fill(event->cent9(), Qy2_west);
    P_Q3_cent_west->Fill(event->cent9(), Q3_west);
    P_Qx3_cent_west->Fill(event->cent9(), Qx3_west);
    P_Qy3_cent_west->Fill(event->cent9(), Qy3_west);
    P_Qx2_cent_RunID_west->Fill(event->runId(), event->cent9(), Qx2_west, 1);
    P_Qy2_cent_RunID_west->Fill(event->runId(), event->cent9(), Qy2_west, 1);
    P_Qx3_cent_RunID_west->Fill(event->runId(), event->cent9(), Qx3_west, 1);
    P_Qy3_cent_RunID_west->Fill(event->runId(), event->cent9(), Qy3_west, 1);
    //----east
    P_Q2_cent_east->Fill(event->cent9(), Q2_east);
    P_Qx2_cent_east->Fill(event->cent9(), Qx2_east);
    P_Qy2_cent_east->Fill(event->cent9(), Qy2_east);
    P_Q3_cent_east->Fill(event->cent9(), Q3_east);
    P_Qx3_cent_east->Fill(event->cent9(), Qx3_east);
    P_Qy3_cent_east->Fill(event->cent9(), Qy3_east);
    P_Qx2_cent_RunID_east->Fill(event->runId(), event->cent9(), Qx2_east, 1);
    P_Qy2_cent_RunID_east->Fill(event->runId(), event->cent9(), Qy2_east, 1);
    P_Qx3_cent_RunID_east->Fill(event->runId(), event->cent9(), Qx3_east, 1);
    P_Qy3_cent_RunID_east->Fill(event->runId(), event->cent9(), Qy3_east, 1);

    //Filling hystograms for Q-Vectors
    H_Q2->Fill(Q2);
    H_Q3->Fill(Q3);
    H_Qx2->Fill(Qx2);
    H_Qy2->Fill(Qy2);
    H_Qx3->Fill(Qy3);
    H_Qy3->Fill(Qy3);
    H_Psi2->Fill(Psi2);
    H_Psi3->Fill(Psi3);
    H_Psi2_west->Fill(Psi2_west);
    H_Psi3_west->Fill(Psi3_west);
    H_Psi2_east->Fill(Psi2_east);
    H_Psi3_east->Fill(Psi3_east);

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  //Saving our files
  TFile *savefile = new TFile(outFileName, "RECREATE");

  H_refMult->Write();
  H_cent9->Write();
  H_VtxZ->Write();
  H_VxVy->Write();
  H_runID->Write();
  H_massSqr->Write();
  H_invBeta->Write();
  H_dEdx->Write();
  H_Pt->Write();
  H_Px->Write();
  H_Py->Write();
  H_Phi->Write();
  H_Eta->Write();
  H_PtEta->Write();
  H_PtPhi->Write();
  H_EtaPhi->Write();
  H_nHits->Write();
  H_dEdxQP->Write();
  H_m2QP->Write();
  H_nSigP_QP->Write();
  H_nSigPi_QP->Write();
  H_nSigK_QP->Write();
  // H_MinvPiTPC->Write();
  // H_MinvKaTPC->Write();
  // H_MinvPPiPosTPC->Write();
  // H_MinvPiTOF->Write();
  // H_MinvKaTOF->Write();
  // H_MinvPPiPosTOF->Write();
  H_Q2->Write();
  H_Q3->Write();
  H_Qx2->Write();
  H_Qy2->Write();
  H_Qx3->Write();
  H_Qy3->Write();
  H_Psi3->Write();
  H_Psi2->Write();
  P_Q2_cent->Write();
  P_Qx2_cent->Write();
  P_Qy2_cent->Write();
  P_Q3_cent->Write();
  P_Qx3_cent->Write();
  P_Qy3_cent->Write();
  P_Qx2_cent_RunID->Write();
  P_Qy2_cent_RunID->Write();
  P_Qx3_cent_RunID->Write();
  P_Qy3_cent_RunID->Write();
  //----west
  H_Psi3_west->Write();
  H_Psi2_west->Write();
  P_Q2_cent_west->Write();
  P_Qx2_cent_west->Write();
  P_Qy2_cent_west->Write();
  P_Q3_cent_west->Write();
  P_Qx3_cent_west->Write();
  P_Qy3_cent_west->Write();
  P_Qx2_cent_RunID_west->Write();
  P_Qy2_cent_RunID_west->Write();
  P_Qx3_cent_RunID_west->Write();
  P_Qy3_cent_RunID_west->Write();
  //----east
  H_Psi3_east->Write();
  H_Psi2_east->Write();
  P_Q2_cent_east->Write();
  P_Qx2_cent_east->Write();
  P_Qy2_cent_east->Write();
  P_Q3_cent_east->Write();
  P_Qx3_cent_east->Write();
  P_Qy3_cent_east->Write();
  P_Qx2_cent_RunID_east->Write();
  P_Qy2_cent_RunID_east->Write();
  P_Qx3_cent_RunID_east->Write();
  P_Qy3_cent_RunID_east->Write();

  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}
