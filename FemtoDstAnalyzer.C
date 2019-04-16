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

  Double_t omega=0,Q2=0,Qx2=0,Qy2=0,Psi2=0;
  Double_t Q3=0,Qx3=0,Qy3=0,Psi3=0;
  const Int_t ngap=4; // количество по eta, размер всех массивов
  Double_t Qx2_east[ngap]={ 0.0 }, Qx2_west[ngap]={ 0.0 };
  Double_t Qx3_east[ngap]={ 0.0 }, Qx3_west[ngap]={ 0.0 };
  Double_t Qy2_east[ngap]={ 0.0 }, Qy2_west[ngap]={ 0.0 };
  Double_t Qy3_east[ngap]={ 0.0 }, Qy3_west[ngap]={ 0.0 };
  Double_t Psi2_west[ngap]={0.0}, Psi2_east[ngap]={0.0};
  Double_t Psi3_west[ngap]={0.0}, Psi3_east[ngap]={0.0};
  Double_t N_east[ngap]={ 0.0 },   N_west[ngap]={ 0.0 };
  Double_t N=0;
  Int_t q=0;
  Double_t n[ngap]={0.0, 0.075, 0.05, 0.5};  
  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0}; 

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

  TH1D *H_Psi2_east[ngap][9];
  TH1D *H_Psi2_west[ngap][9];
  TH1D *H_Psi3_east[ngap][9];
  TH1D *H_Psi3_west[ngap][9];

  for (Int_t i=0; i < ngap; i++){
    for(Int_t j=0; j<9; j++){
      H_Psi2_west[i][j] = new TH1D (Form("H_Psi2_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 2 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{2},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi2_east[i][j] = new TH1D (Form("H_Psi2_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 2 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{2},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_west[i][j] = new TH1D (Form("H_Psi3_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_east[i][j] = new TH1D (Form("H_Psi3_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
    }
  }
      
   //Q-vectors WEST hystograms:
  // TH1D *H_Qx2_west = new TH1D("H_Qx2_west","Q_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qy2_west = new TH1D("H_Qy2_west","Q_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Q2_west = new TH1D("H_Q2_west","Q_{2}-vector; Q_{2},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qx3_west = new TH1D("H_Qx3_west","Q_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qy3_west = new TH1D("H_Qy3_west","Q_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q3_west = new TH1D("H_Q3_west","Q_{3}-vector; Q_{3},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Psi2_west = new TH1D("H_Psi2_west","#psi of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  // TH1D *H_Psi3_west = new TH1D("H_Psi3_west","#psi of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  //Q-vectors EAST hystograms:
  // TH1D *H_Qx2_east = new TH1D("H_Qx2_east","Q_{2}_{x} projection; Q_{2}_{x},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qy2_east = new TH1D("H_Qy2_east","Q_{2}_{y} projection; Q_{2}_{y},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Q2_east = new TH1D("H_Q2_east","Q_{2}-vector; Q_{2},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qx3_east = new TH1D("H_Qx3_east","Q_{3}_{x} projection; Q_{3}_{x},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Qy3_east = new TH1D("H_Qy3_east","Q_{3}_{y} projection; Q_{3}_{y},[GeV/c];N_{count}",500,-60, 60);
  TH1D *H_Q3_east = new TH1D("H_Q3_east","Q_{3}-vector; Q_{3},[GeV/c];N_{count}",500,-60, 60);
  // TH1D *H_Psi2_east = new TH1D("H_Psi2_east","#psi of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  // TH1D *H_Psi3_east = new TH1D("H_Psi3_east","#psi of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );

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
  TProfile2D *P_Qx2_cent_RunID = new TProfile2D("P_Qx2_cent_RunID","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy2_cent_RunID = new TProfile2D("P_Qy2_cent_RunID","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qx3_cent_RunID = new TProfile2D("P_Qx3_cent_RunID","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  TProfile2D *P_Qy3_cent_RunID = new TProfile2D("P_Qy3_cent_RunID","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // //-----west
  // TProfile2D *P_Qx2_cent_RunID_west = new TProfile2D("P_Qx2_cent_RunID_west","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qy2_cent_RunID_west = new TProfile2D("P_Qy2_cent_RunID_west","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qx3_cent_RunID_west = new TProfile2D("P_Qx3_cent_RunID_west","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qy3_cent_RunID_west = new TProfile2D("P_Qy3_cent_RunID_west","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );

  TProfile2D *P_Qx2_cent_RunID_west[ngap];
  TProfile2D *P_Qy2_cent_RunID_west[ngap];
  TProfile2D *P_Qx3_cent_RunID_west[ngap];
  TProfile2D *P_Qy3_cent_RunID_west[ngap];
  
  for (int i=0; i<ngap; i++){
    P_Qx2_cent_RunID_west[i] = new TProfile2D (Form("P_Qx2_cent_RunID_west_n=%.2f",2.0*n[i]), Form("Profile of Q_{2}_{x} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qy2_cent_RunID_west[i] = new TProfile2D (Form("P_Qy2_cent_RunID_west_n=%.2f",2.0*n[i]), Form("Profile of Q_{2}_{y} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qx3_cent_RunID_west[i] = new TProfile2D (Form("P_Qx3_cent_RunID_west_n=%.2f",2.0*n[i]), Form("Profile of Q_{3}_{x} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qy3_cent_RunID_west[i] = new TProfile2D (Form("P_Qy3_cent_RunID_west_n=%.2f",2.0*n[i]), Form("Profile of Q_{3}_{y} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }

  //-----east
  // TProfile2D *P_Qx2_cent_RunID_east = new TProfile2D("P_Qx2_cent_RunID_east","Profile of Q_{2}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qy2_cent_RunID_east = new TProfile2D("P_Qy2_cent_RunID_east","Profile of Q_{2}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qx3_cent_RunID_east = new TProfile2D("P_Qx3_cent_RunID_east","Profile of Q_{3}_{x} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );
  // TProfile2D *P_Qy3_cent_RunID_east = new TProfile2D("P_Qy3_cent_RunID_east","Profile of Q_{3}_{y} on Centrality vs RunID",0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5 );

  TProfile2D *P_Qx2_cent_RunID_east[ngap];
  TProfile2D *P_Qy2_cent_RunID_east[ngap];
  TProfile2D *P_Qx3_cent_RunID_east[ngap];
  TProfile2D *P_Qy3_cent_RunID_east[ngap];
  
  for (int i=0; i<ngap; i++){
    P_Qx2_cent_RunID_east[i] = new TProfile2D (Form("P_Qx2_cent_RunID_east_n=%.2f",2.0*n[i]), Form("Profile of Q_{2}_{x} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qy2_cent_RunID_east[i] = new TProfile2D (Form("P_Qy2_cent_RunID_east_n=%.2f",2.0*n[i]), Form("Profile of Q_{2}_{y} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qx3_cent_RunID_east[i] = new TProfile2D (Form("P_Qx3_cent_RunID_east_n=%.2f",2.0*n[i]), Form("Profile of Q_{3}_{x} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
    P_Qy3_cent_RunID_east[i] = new TProfile2D (Form("P_Qy3_cent_RunID_east_n=%.2f",2.0*n[i]), Form("Profile of Q_{3}_{y} on Centrality vs RunID, #eta=%.2f;RunId; Centrality",2.0*n[i]),0.01e6, 12170000.0-0.0005, 12180000.0+0.0005, 9, -0.5, 8.5);
  }
  

   


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
        if ((iEvent+1) % 10000 == 0){
      std::cout << "Working on event #[" << (iEvent + 1)
                << "/" << events2read << "]" << std::endl;
    }

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
    if( TMath::Abs(pow( pVtx(0), 2 ) + pow( pVtx(1) , 2)) > 2. ) continue;
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
    for (int i = 0; i < ngap; i++)
    {
      Qx2_east[i]= 0.0; Qx2_west[i]=0.0;
      Qx3_east[i]=0.0; Qx3_west[i]=0.0;
      Qy2_east[i]=0.0; Qy2_west[i]=0.0;
      Qy3_east[i]=0.0; Qy3_west[i]=0.0;
      Psi2_west[i]=0.0; Psi2_east[i]=0.0;
      Psi3_west[i]=0.0; Psi3_east[i]=0.0;
      N_east[i]= 0; N_west[i]=0;
    }
    N=0;

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
      if( femtoTrack->gMom().Mag() < 0.1 || femtoTrack->gDCA(pVtx).Mag() > 1. ) {
         continue;
      }
      /*//Для индитификации частиц
      if( TMath::Abs( femtoTrack -> eta() ) > 1.0 || femtoTrack -> nHits() < 15 || 
          femtoTrack -> p() < 0.15 || femtoTrack -> p() > 5.0) {
         continue;
      }
      */
      // для заряженных адронов
      if( TMath::Abs( femtoTrack -> eta() ) > 1.0 ||
          femtoTrack -> nHits() < 15 ||
          femtoTrack -> pt() < 0.2 || 
          femtoTrack -> pt() > 2.0 ||
          femtoTrack -> p() < 0.15 || 
          femtoTrack -> p() > 5.0) {
         continue;
      }
      
      

      omega=femtoTrack->pt();
      //Qx and Qy
      Qx2+= omega*cos(2*(femtoTrack->phi()) );
      Qy2+= omega*sin(2*(femtoTrack->phi()) );
      Qx3+= omega*cos(3*(femtoTrack->phi()) );
      Qy3+= omega*sin(3*(femtoTrack->phi()) );
      N++;
      for(Int_t i=0; i<ngap; i++){
        if(femtoTrack->eta()>n[i]){
          Qx2_east[i]+= omega*cos(2*(femtoTrack->phi()) );
          Qy2_east[i]+= omega*sin(2*(femtoTrack->phi()) );
          Qx3_east[i]+= omega*cos(3*(femtoTrack->phi()) );
          Qy3_east[i]+= omega*sin(3*(femtoTrack->phi()) );
          N_east[i]++;
        }

        if(femtoTrack->eta()< -n[i]){
          Qx2_west[i]+= omega*cos(2*(femtoTrack->phi()) );
          Qy2_west[i]+= omega*sin(2*(femtoTrack->phi()) );
          Qx3_west[i]+= omega*cos(3*(femtoTrack->phi()) );
          Qy3_west[i]+= omega*sin(3*(femtoTrack->phi()) );
          N_west[i]++;
        }
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


    if( N != 0){
        Qx2= Qx2 / N;
        Qy2= Qy2 / N;
        Qx3= Qx3 / N;
        Qy3= Qy3 / N;
    }
    Psi2=1.0/2.0*(TMath::ATan2(Qy2, Qx2)); //(TMath::Sqrt(event->refMult()))
    Psi3=1.0/3.0*(TMath::ATan2(Qy3, Qx3));
    P_Qx2_cent_RunID->Fill(event->runId(), event->cent9(), Qx2, 1);
    P_Qy2_cent_RunID->Fill(event->runId(), event->cent9(), Qy2, 1);
    P_Qx3_cent_RunID->Fill(event->runId(), event->cent9(), Qx3, 1);
    P_Qy3_cent_RunID->Fill(event->runId(), event->cent9(), Qy3, 1);

    //Filling hystograms for Q-Vectors
    H_Q2->Fill(Q2);
    H_Q3->Fill(Q3);
    H_Qx2->Fill(Qx2);
    H_Qy2->Fill(Qy2);
    H_Qx3->Fill(Qy3);
    H_Qy3->Fill(Qy3);
    H_Psi2->Fill(Psi2);
    H_Psi3->Fill(Psi3);

    for(Int_t i=0; i<ngap; i++){
      if( N_west[i] != 0 ){
        Qx2_west[i] = Qx2_west[i] / N_west[i]; 
        Qy2_west[i] = Qy2_west[i] / N_west[i];
        Qx3_west[i] = Qx3_west[i] / N_west[i];
        Qy3_west[i] = Qy3_west[i] / N_west[i];
      }
      if( N_east[i] != 0 ){
        Qx2_east[i] = Qx2_east[i] / N_east[i]; 
        Qy2_east[i] = Qy2_east[i] / N_east[i];
        Qx3_east[i] = Qx3_east[i] / N_east[i];
        Qy3_east[i] = Qy3_east[i] / N_east[i];
      }
      //------west arm Q
      if( N_west[i] != 0 ){
        Psi2_west[i]=1.0/2.0*(TMath::ATan2(Qy2_west[i], Qx2_west[i])); //(TMath::Sqrt(event->refMult()))
        Psi3_west[i]=1.0/3.0*(TMath::ATan2(Qy3_west[i], Qx3_west[i]));
        //Q2_west=TMath::Sqrt(pow(Qx2_west[i],2)+pow(Qy2_west[i],2));
        //Q3_west=TMath::Sqrt(pow(Qx3_west[i],2)+pow(Qy3_west[i],2));
      }
        //-----east Q
      if( N_east[i] != 0 ){
        Psi2_east[i]=1.0/2.0*(TMath::ATan2(Qy2_east[i], Qx2_east[i])); //(TMath::Sqrt(event->refMult()))
        Psi3_east[i]=1.0/3.0*(TMath::ATan2(Qy3_east[i], Qx3_east[i])); 
        //Q2_east=TMath::Sqrt(pow(Qx2_east[i],2)+pow(Qy2_east[i],2));
        //Q3_east=TMath::Sqrt(pow(Qx3_east[i],2)+pow(Qy3_east[i],2));
      }
      //Fillig Profiles
      //----west
      P_Qx2_cent_RunID_west[i]->Fill(event->runId(), event->cent9(), Qx2_west[i], 1);
      P_Qy2_cent_RunID_west[i]->Fill(event->runId(), event->cent9(), Qy2_west[i], 1);
      P_Qx3_cent_RunID_west[i]->Fill(event->runId(), event->cent9(), Qx3_west[i], 1);
      P_Qy3_cent_RunID_west[i]->Fill(event->runId(), event->cent9(), Qy3_west[i], 1);
      //----east
      P_Qx2_cent_RunID_east[i]->Fill(event->runId(), event->cent9(), Qx2_east[i], 1);
      P_Qy2_cent_RunID_east[i]->Fill(event->runId(), event->cent9(), Qy2_east[i], 1);
      P_Qx3_cent_RunID_east[i]->Fill(event->runId(), event->cent9(), Qx3_east[i], 1);
      P_Qy3_cent_RunID_east[i]->Fill(event->runId(), event->cent9(), Qy3_east[i], 1);
      
        H_Psi2_west[i][event->cent9()]->Fill(Psi2_west[i]);
        H_Psi3_west[i][event->cent9()]->Fill(Psi3_west[i]);
        H_Psi2_east[i][event->cent9()]->Fill(Psi2_east[i]);
        H_Psi3_east[i][event->cent9()]->Fill(Psi3_east[i]);
      
    }

    
      

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
  P_Qx2_cent_RunID->Write();
  P_Qy2_cent_RunID->Write();
  P_Qx3_cent_RunID->Write();
  P_Qy3_cent_RunID->Write();
  for(Int_t i=0; i<ngap; i++){
    //----west
    for(Int_t j=0; j<9; j++){
      H_Psi3_west[i][j]->Write();
      H_Psi2_west[i][j]->Write();
    }
    P_Qx2_cent_RunID_west[i]->Write();
    P_Qy2_cent_RunID_west[i]->Write();
    P_Qx3_cent_RunID_west[i]->Write();
    P_Qy3_cent_RunID_west[i]->Write();
    //----east
    for(Int_t j=0; j<9; j++){
      H_Psi3_east[i][j]->Write();
      H_Psi2_east[i][j]->Write();
    }
    P_Qx2_cent_RunID_east[i]->Write();
    P_Qy2_cent_RunID_east[i]->Write();
    P_Qx3_cent_RunID_east[i]->Write();
    P_Qy3_cent_RunID_east[i]->Write();
  }

  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}
