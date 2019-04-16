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
void FemtoDstAnalyzerFlat(const Char_t *inFile = "AuAu27GeV/AuAu27_ar.list", const Char_t *outFileName = "outFile27.root",const Char_t *flatFileName="Result27Recenter.root") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  //Saving our files
  TFile *savefile = new TFile(outFileName, "RECREATE");

  Double_t omega=0, Qx2=0, Qy2=0, Qx2_recenter=0, Qy2_recenter=0, Psi2=0, Psi2_recenter=0, delta_Psi2=0, Psi2_flat=0;
  Double_t Qx3=0, Qy3=0, Qx3_recenter=0, Qy3_recenter=0, Psi3=0, Psi3_recenter=0, delta_Psi3=0, Psi3_flat=0;
  const Int_t ngap=4; // количество по eta, размер всех массивов
  Double_t Qx2_east[ngap]={ 0.0 }, Qx2_west[ngap]={ 0.0 };
  Double_t Qx3_east[ngap]={ 0.0 }, Qx3_west[ngap]={ 0.0 };
  Double_t Qy2_east[ngap]={ 0.0 }, Qy2_west[ngap]={ 0.0 };
  Double_t Qy3_east[ngap]={ 0.0 }, Qy3_west[ngap]={ 0.0 };
  Double_t Qx2_recenter_east[ngap]={ 0.0 }, Qx2_recenter_west[ngap]={ 0.0 };
  Double_t Qy2_recenter_east[ngap]={ 0.0 }, Qy2_recenter_west[ngap]={ 0.0 };
  Double_t Qx3_recenter_east[ngap]={ 0.0 }, Qx3_recenter_west[ngap]={ 0.0 };
  Double_t Qy3_recenter_east[ngap]={ 0.0 }, Qy3_recenter_west[ngap]={ 0.0 };
  Double_t Psi2_west[ngap]={0.0}, Psi2_east[ngap]={0.0};
  Double_t Psi3_west[ngap]={0.0}, Psi3_east[ngap]={0.0};
  Double_t Psi2_recenter_west[ngap]={0.0}, Psi2_recenter_east[ngap]={0.0};
  Double_t Psi3_recenter_west[ngap]={0.0}, Psi3_recenter_east[ngap]={0.0}; 

  Double_t Psi2_flat_west[ngap]={0.0}, Psi3_flat_west[ngap]={0.0}, delta_Psi2_west[ngap]={0.0}, delta_Psi3_west[ngap]={0.0};
  Double_t Psi2_flat_east[ngap]={0.0}, Psi3_flat_east[ngap]={0.0}, delta_Psi2_east[ngap]={0.0}, delta_Psi3_east[ngap]={0.0};
  Double_t N_east[ngap]={ 0.0 },   N_west[ngap]={ 0.0 };
  Double_t N=0;
  Int_t q=0;
  Double_t n[ngap]={0.0, 0.075, 0.05, 0.5};  
  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0}; 

  TH2D *H_dEdxQP = new TH2D("dEdxQP","dE/dx vs q*p; Q*P,[GeV/c]; dE/dx,[a.u.]",1000,-1.5,1.5, 1000,0.0,14.0);
  TH2D *H_m2qpt = new TH2D("H_m2qpt","m^{2} vs q*P_{t}; Q*P_{t},[GeV/c]; m^{2},[(GeV/c)^{2}]", 400,-2.1,2.1, 400,-0.5,1.5);
  TH2D *H_nSigmaPion = new TH2D("nSigmaPionQP","nSigmaPion vs q*p;Q*P;nSigmaPion",200,-3,3, 200,-10,10);
  TH2D *H_nSigmaKaon = new TH2D("nSigmaKaonQP","nSigmaKaon vs q*p;Q*P;nSigmaKaon",200,-3,3, 200,-20,25);
  TH2D *H_nSigmaProton = new TH2D("nSigmaProtonQP","nSigmaProton vs q*p;Q*P;nSigmaProton",200,-3,3, 200,-25,25);

  TH1D *H_Psi2 = new TH1D("H_Psi2","#psi_{2} of Event plane of 2 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3 = new TH1D("H_Psi3","#psi_{3} of Event plane of 3 harmonic; #psi,[rad]",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi2_recenter = new TH1D("H_Psi2_recenter","#psi^{recenter}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_recenter = new TH1D("H_Psi3_recenter","#psi^{recenter}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  TH1D *H_Psi2_flat = new TH1D("H_Psi2_flat","#psi^{rec+flat}_{2} of Event plane of 2 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/2.0-1, TMath::Pi()/2.0+1 );
  TH1D *H_Psi3_flat = new TH1D("H_Psi3_flat","#psi^{rec+flat}_{3} of Event plane of 3 harmonic; #psi,[rad]; N_{count}",500,-TMath::Pi()/3.0-1, TMath::Pi()/3.0+1 );
  

  TH1D *H_Psi2_east[ngap][9];
  TH1D *H_Psi2_west[ngap][9];
  TH1D *H_Psi3_east[ngap][9];
  TH1D *H_Psi3_west[ngap][9];
  TH1D *H_Psi2_recenter_west[ngap][9];
  TH1D *H_Psi2_recenter_east[ngap][9];
  TH1D *H_Psi3_recenter_west[ngap][9];
  TH1D *H_Psi3_recenter_east[ngap][9];
  TH1D *H_Psi2_flat_west[ngap][9];
  TH1D *H_Psi2_flat_east[ngap][9];
  TH1D *H_Psi3_flat_west[ngap][9];
  TH1D *H_Psi3_flat_east[ngap][9];
  TProfile *P_square_resolution2[ngap];
  TProfile *P_square_resolution3[ngap];

  for (Int_t i=0; i < ngap; i++){
    for(Int_t j=0; j<9; j++){
      H_Psi2_west[i][j] = new TH1D (Form("H_Psi2_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 2 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{2},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi2_east[i][j] = new TH1D (Form("H_Psi2_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 2 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{2},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_west[i][j] = new TH1D (Form("H_Psi3_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_east[i][j] = new TH1D (Form("H_Psi3_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);

      H_Psi2_recenter_west[i][j] = new TH1D (Form("H_Psi2_recenter_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi2_recenter_east[i][j] = new TH1D (Form("H_Psi2_recenter_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_recenter_west[i][j] = new TH1D (Form("H_Psi3_recenter_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_recenter_east[i][j] = new TH1D (Form("H_Psi3_recenter_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);

      H_Psi2_flat_west[i][j] = new TH1D (Form("H_Psi2_flat_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec+flat} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi2_flat_east[i][j] = new TH1D (Form("H_Psi2_flat_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec+flat} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_flat_west[i][j] = new TH1D (Form("H_Psi3_flat_west_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec+flat} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
      H_Psi3_flat_east[i][j] = new TH1D (Form("H_Psi3_flat_east_n=%.2f_centr9=%i",2.0*n[i],cent[j]), Form("#psi^{rec+flat} of Event plane of 3 harmonic(west) #eta-gap=%.2f, centr %i - %i %% ;#Psi_{3},[rad];Counts",2.0*n[i],cent[j+1],cent[j]),500,-TMath::Pi()/2.0-0.1,TMath::Pi()/2.0+0.1);
    }

    //Profiles for resolution stage
    P_square_resolution2[i] = new TProfile(Form("P_square_resolution2_n=%i",i+1),Form("Profile of <cos(2*(#psi_{2,#pm}-#psi_{2,#mp}))>, #eta-gap=%.2f; Centrality bins",2.0*n[i]), 9, -0.5, 8.5);
    P_square_resolution3[i] = new TProfile(Form("P_square_resolution3_n=%i",i+1),Form("Profile of <cos(2*(#psi_{3,#pm}-#psi_{3,#mp}))>, #eta-gap=%.2f; Centrality bins",2.0*n[i]), 9, -0.5, 8.5);
  }
  
  
  //---------------Чтение файла
  TFile *readf = new TFile(flatFileName,"READ");
  TProfile2D *P3_Qx2 = (TProfile2D*)readf->Get("P2_Qx2_cent_RunID");
  TProfile2D *P3_Qy2 = (TProfile2D*)readf->Get("P2_Qy2_cent_RunID");
  TProfile2D *P3_Qx3 = (TProfile2D*)readf->Get("P2_Qx3_cent_RunID");
  TProfile2D *P3_Qy3 = (TProfile2D*)readf->Get("P2_Qy3_cent_RunID");
  // Profile for cos and sin at v2
  TProfile2D *p_sin_v2[4];
  TProfile2D *p_cos_v2[4];
  TProfile2D *p_sin_v3[4];
  TProfile2D *p_cos_v3[4];
  for (int i=0; i<4; i++){
    p_sin_v2[i] = (TProfile2D*)readf -> Get(Form("sin_v2_prof_k=%i",i+1));
    p_cos_v2[i] = (TProfile2D*)readf -> Get(Form("cos_v2_prof_k=%i",i+1));
    p_sin_v3[i] = (TProfile2D*)readf -> Get(Form("sin_v3_prof_k=%i",i+1));
    p_cos_v3[i] = (TProfile2D*)readf -> Get(Form("cos_v3_prof_k=%i",i+1));
  }  



  TProfile2D *P3_Qx2_west[ngap];
  TProfile2D *P3_Qy2_west[ngap];
  TProfile2D *P3_Qx3_west[ngap];
  TProfile2D *P3_Qy3_west[ngap];
  TProfile2D *P3_Qx2_east[ngap];
  TProfile2D *P3_Qy2_east[ngap];
  TProfile2D *P3_Qx3_east[ngap];
  TProfile2D *P3_Qy3_east[ngap];
  //Профайлы синусов и косинусов с учётом гэпов
  TProfile2D *p_sin_v2_west[4][ngap];
  TProfile2D *p_cos_v2_west[4][ngap];
  TProfile2D *p_sin_v3_west[4][ngap];
  TProfile2D *p_cos_v3_west[4][ngap];

  TProfile2D *p_sin_v2_east[4][ngap];
  TProfile2D *p_cos_v2_east[4][ngap];
  TProfile2D *p_sin_v3_east[4][ngap];
  TProfile2D *p_cos_v3_east[4][ngap];


  for (Int_t i=0; i < ngap; i++){
    //_west
    P3_Qx2_west[i] = (TProfile2D*)readf->Get(Form("P2_Qx2_cent_RunID_west_n=%.2f",2.0*n[i]));
    P3_Qy2_west[i] = (TProfile2D*)readf->Get(Form("P2_Qy2_cent_RunID_west_n=%.2f",2.0*n[i]));
    P3_Qx3_west[i] = (TProfile2D*)readf->Get(Form("P2_Qx3_cent_RunID_west_n=%.2f",2.0*n[i]));
    P3_Qy3_west[i] = (TProfile2D*)readf->Get(Form("P2_Qy3_cent_RunID_west_n=%.2f",2.0*n[i]));
    //_east
    P3_Qx2_east[i] = (TProfile2D*)readf->Get(Form("P2_Qx2_cent_RunID_east_n=%.2f",2.0*n[i]));
    P3_Qy2_east[i] = (TProfile2D*)readf->Get(Form("P2_Qy2_cent_RunID_east_n=%.2f",2.0*n[i]));
    P3_Qx3_east[i] = (TProfile2D*)readf->Get(Form("P2_Qx3_cent_RunID_east_n=%.2f",2.0*n[i]));
    P3_Qy3_east[i] = (TProfile2D*)readf->Get(Form("P2_Qy3_cent_RunID_east_n=%.2f",2.0*n[i]));
    //WEST
    for (int k=0; k<4; k++){
      p_sin_v2_west[k][i] = (TProfile2D*)readf -> Get(Form("sin_v2_west_prof_k=%i_n=%i",k+1,i+1));
      p_cos_v2_west[k][i] = (TProfile2D*)readf -> Get(Form("cos_v2_west_prof_k=%i_n=%i",k+1,i+1));
      p_sin_v3_west[k][i] = (TProfile2D*)readf -> Get(Form("sin_v3_west_prof_k=%i_n=%i",k+1,i+1));
      p_cos_v3_west[k][i] = (TProfile2D*)readf -> Get(Form("cos_v3_west_prof_k=%i_n=%i",k+1,i+1));
    //EAST
      p_sin_v2_east[k][i] = (TProfile2D*)readf -> Get(Form("sin_v2_east_prof_k=%i_n=%i",k+1,i+1));
      p_cos_v2_east[k][i] = (TProfile2D*)readf -> Get(Form("cos_v2_east_prof_k=%i_n=%i",k+1,i+1));
      p_sin_v3_east[k][i] = (TProfile2D*)readf -> Get(Form("sin_v3_east_prof_k=%i_n=%i",k+1,i+1));
      p_cos_v3_east[k][i] = (TProfile2D*)readf -> Get(Form("cos_v3_east_prof_k=%i_n=%i",k+1,i+1));
    }
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


  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

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

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
 
    //Q-vector cleaning in new event
    Qx2=0; Qy2=0; Psi2=0; Qx3=0; Qy3=0; Psi3=0;
    Qx2_recenter=0; Qy2_recenter=0; Psi2_recenter=0;
    Qx3_recenter=0; Qy3_recenter=0; Psi3_recenter=0;
    delta_Psi3=0; Psi3_flat=0;
    delta_Psi2=0; Psi2_flat=0;
    for (int i = 0; i < ngap; i++)
    {
      Qx2_east[i]= 0.0; Qx2_west[i]=0.0;
      Qx3_east[i]=0.0; Qx3_west[i]=0.0;
      Qy2_east[i]=0.0; Qy2_west[i]=0.0;
      Qy3_east[i]=0.0; Qy3_west[i]=0.0;
      Qx2_recenter_east[i]=0.0; Qx2_recenter_west[i]=0.0;
      Qy2_recenter_east[i]=0.0; Qy2_recenter_west[i]=0.0;
      Qx3_recenter_east[i]=0.0; Qx3_recenter_west[i]=0.0;
      Qy3_recenter_east[i]=0.0; Qy3_recenter_west[i]=0.0;
      Psi2_west[i]=0.0; Psi2_east[i]=0.0;
      Psi3_west[i]=0.0; Psi3_east[i]=0.0;
      Psi2_recenter_west[i]=0.0; Psi2_recenter_east[i]=0.0;
      Psi3_recenter_west[i]=0.0; Psi3_recenter_east[i]=0.0;
      Psi2_flat_west[i]=0.0; Psi2_flat_east[i]=0.0;
      Psi3_flat_west[i]=0.0; Psi3_flat_east[i]=0.0;
      delta_Psi2_west[i]=0; delta_Psi3_west[i]=0;
      delta_Psi2_east[i]=0; delta_Psi3_east[i]=0;
      N_east[i]= 0; N_west[i]=0;
    }
    N=0;

    // Track loop
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

      H_dEdxQP->Fill(femtoTrack->p()*femtoTrack->charge(),femtoTrack->dEdx()*1e6);
      H_nSigmaPion->Fill(femtoTrack->p()*femtoTrack->charge(),femtoTrack->nSigmaPion());
      H_nSigmaKaon->Fill(femtoTrack->p()*femtoTrack->charge(),femtoTrack->nSigmaKaon());
      H_nSigmaProton->Fill(femtoTrack->p()*femtoTrack->charge(),femtoTrack->nSigmaProton());

      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {
        H_m2qpt->Fill(femtoTrack->pt()*femtoTrack->charge(),femtoTrack->massSqr());

      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    if( N != 0){
        Qx2= Qx2 / N;
        Qy2= Qy2 / N;
        Qx3= Qx3 / N;
        Qy3= Qy3 / N;
    }
    Psi2=1.0/2.0*(TMath::ATan2(Qy2, Qx2));
    Psi3=1.0/3.0*(TMath::ATan2(Qy3, Qx3));
    Qx2_recenter=Qx2- P3_Qx2->GetBinContent(P3_Qx2->FindBin(event->runId(), event->cent9()));
    Qy2_recenter=Qy2- P3_Qy2->GetBinContent(P3_Qy2->FindBin(event->runId(), event->cent9()));
    Qx3_recenter=Qx3- P3_Qx3->GetBinContent(P3_Qx3->FindBin(event->runId(), event->cent9()));
    Qy3_recenter=Qy3- P3_Qy3->GetBinContent(P3_Qy3->FindBin(event->runId(), event->cent9()));
    Psi2_recenter=1.0/2.0*(TMath::ATan2(Qy2_recenter, Qx2_recenter)); 
    Psi3_recenter=1.0/3.0*(TMath::ATan2(Qy3_recenter, Qx3_recenter));
    for (int i=0; i<4; i++){
      delta_Psi2 += 2.0/( (Double_t)i+1.0)*(p_cos_v2[i]->GetBinContent(p_cos_v2[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*2*Psi2_recenter)
      -p_sin_v2[i]->GetBinContent(p_sin_v2[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*2*Psi2_recenter) );
    }
    for (int i=0; i<4; i++){
      delta_Psi3 += 2.0/( (Double_t)i+1.0)*(p_cos_v3[i]->GetBinContent(p_cos_v3[i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*3*Psi3_recenter)
      -p_sin_v3[i]->GetBinContent(p_sin_v3[i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*3*Psi3_recenter) );
    }
    Psi2_flat = Psi2_recenter + (1.0/2.0)*delta_Psi2;
    Psi3_flat = Psi3_recenter + (1.0/3.0)*delta_Psi3;

    //Filling hystograms for Q-Vectors
    H_Psi2->Fill(Psi2);
    H_Psi3->Fill(Psi3);
    H_Psi2_recenter->Fill(Psi2_recenter);
    H_Psi3_recenter->Fill(Psi3_recenter);
    H_Psi2_flat->Fill(Psi2_flat);
    H_Psi3_flat->Fill(Psi3_flat);





    //--------------------Заполнение для west и east----------------
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
          Psi2_west[i]=1.0/2.0*(TMath::ATan2(Qy2_west[i], Qx2_west[i]));
          Psi3_west[i]=1.0/3.0*(TMath::ATan2(Qy3_west[i], Qx3_west[i]));
      }
      //-----east Q
      if( N_east[i] != 0 ){
          Psi2_east[i]=1.0/2.0*(TMath::ATan2(Qy2_east[i], Qx2_east[i]));
          Psi3_east[i]=1.0/3.0*(TMath::ATan2(Qy3_east[i], Qx3_east[i]));
      }

      //----_west
      if( N_west[i] != 0 ){
        Qx2_recenter_west[i]=Qx2_west[i] - P3_Qx2_west[i]->GetBinContent(P3_Qx2_west[i]->FindBin(event->runId(), event->cent9()));
        Qy2_recenter_west[i]=Qy2_west[i] - P3_Qy2_west[i]->GetBinContent(P3_Qy2_west[i]->FindBin(event->runId(), event->cent9()));
        Qx3_recenter_west[i]=Qx3_west[i] - P3_Qx3_west[i]->GetBinContent(P3_Qx3_west[i]->FindBin(event->runId(), event->cent9()));
        Qy3_recenter_west[i]=Qy3_west[i] - P3_Qy3_west[i]->GetBinContent(P3_Qy3_west[i]->FindBin(event->runId(), event->cent9()));
      }
      //----east
      if( N_east[i] != 0 ){
        Qx2_recenter_east[i]=Qx2_east[i] - P3_Qx2_east[i]->GetBinContent(P3_Qx2_east[i]->FindBin(event->runId(), event->cent9()));
        Qy2_recenter_east[i]=Qy2_east[i] - P3_Qy2_east[i]->GetBinContent(P3_Qy2_east[i]->FindBin(event->runId(), event->cent9()));
        Qx3_recenter_east[i]=Qx3_east[i] - P3_Qx3_east[i]->GetBinContent(P3_Qx3_east[i]->FindBin(event->runId(), event->cent9()));
        Qy3_recenter_east[i]=Qy3_east[i] - P3_Qy3_east[i]->GetBinContent(P3_Qy3_east[i]->FindBin(event->runId(), event->cent9()));
      }

      //_west
      if( N_west[i] != 0 ){
        Psi2_recenter_west[i]=1.0/2.0*(TMath::ATan2(Qy2_recenter_west[i], Qx2_recenter_west[i])); 
        Psi3_recenter_west[i]=1.0/3.0*(TMath::ATan2(Qy3_recenter_west[i], Qx3_recenter_west[i]));
      }
      //_east
      if( N_east[i] != 0 ){
        Psi2_recenter_east[i]=1.0/2.0*(TMath::ATan2(Qy2_recenter_east[i], Qx2_recenter_east[i])); 
        Psi3_recenter_east[i]=1.0/3.0*(TMath::ATan2(Qy3_recenter_east[i], Qx3_recenter_east[i]));
      }

      //WEST
      if( N_west[i] != 0 ){
        for (int k=0; k<4; k++){
          delta_Psi2_west[i] += 2.0/( (Double_t)i+1.0)*(p_cos_v2_west[k][i]->GetBinContent(p_cos_v2_west[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*2*Psi2_recenter_west[i])
          -p_sin_v2_west[k][i]->GetBinContent(p_sin_v2_west[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*2*Psi2_recenter_west[i]) );

          delta_Psi3_west[i] += 2.0/( (Double_t)i+1.0)*(p_cos_v3_west[k][i]->GetBinContent(p_cos_v3_west[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*3*Psi3_recenter_west[i])
          -p_sin_v3_west[k][i]->GetBinContent(p_sin_v3_west[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*3*Psi3_recenter_west[i]) );
        }
          Psi2_flat_west[i] = Psi2_recenter_west[i] + (1.0/2.0)*delta_Psi2_west[i];
          Psi3_flat_west[i] = Psi3_recenter_west[i] + (1.0/3.0)*delta_Psi3_west[i];
      }

      //EAST
      if( N_east[i] != 0 ){
        for (int k=0; k<4; k++){
          delta_Psi2_east[i] += 2.0/( (Double_t)i+1.0)*(p_cos_v2_east[k][i]->GetBinContent(p_cos_v2_east[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*2*Psi2_recenter_east[i])
          -p_sin_v2_east[k][i]->GetBinContent(p_sin_v2_east[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*2*Psi2_recenter_east[i]) );

          delta_Psi3_east[i] += 2.0/( (Double_t)i+1.0)*(p_cos_v3_east[k][i]->GetBinContent(p_cos_v3_east[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Sin((i+1)*3*Psi3_recenter_east[i])
          -p_sin_v3_east[k][i]->GetBinContent(p_sin_v3_east[k][i]->FindBin(event->runId(), event->cent9()))*TMath::Cos((i+1)*3*Psi3_recenter_east[i]) );
        }
          Psi2_flat_east[i] = Psi2_recenter_east[i] + (1.0/2.0)*delta_Psi2_east[i];
          Psi3_flat_east[i] = Psi3_recenter_east[i] + (1.0/3.0)*delta_Psi3_east[i];
      }

      //_west
      if( N_west[i] != 0 ){
        H_Psi2_west[i][event->cent9()]->Fill(Psi2_west[i]);
        H_Psi3_west[i][event->cent9()]->Fill(Psi3_west[i]);
        H_Psi2_recenter_west[i][event->cent9()]->Fill(Psi2_recenter_west[i]);
        H_Psi3_recenter_west[i][event->cent9()]->Fill(Psi3_recenter_west[i]);
        H_Psi2_flat_west[i][event->cent9()]->Fill(Psi2_flat_west[i]);
        H_Psi3_flat_west[i][event->cent9()]->Fill(Psi3_flat_west[i]);
      }
      //_east
      if( N_east[i] != 0 ){
        H_Psi2_east[i][event->cent9()]->Fill(Psi2_east[i]);
        H_Psi3_east[i][event->cent9()]->Fill(Psi3_east[i]);
        H_Psi2_recenter_east[i][event->cent9()]->Fill(Psi2_recenter_east[i]);
        H_Psi3_recenter_east[i][event->cent9()]->Fill(Psi3_recenter_east[i]);
        H_Psi2_flat_east[i][event->cent9()]->Fill(Psi2_flat_east[i]);
        H_Psi3_flat_east[i][event->cent9()]->Fill(Psi3_flat_east[i]);
      }
      if(N_west[i]!=0 && N_east[i]!=0){
        P_square_resolution2[i]->Fill(event->cent9(), TMath::Cos(2*( Psi2_flat_west[i] - Psi2_flat_east[i] )) );
        P_square_resolution3[i]->Fill(event->cent9(), TMath::Cos(3*( Psi3_flat_west[i] - Psi3_flat_east[i] )) );
      }
    }

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)




  // H_Psi2->Write();
  // H_Psi3->Write();
  // H_Psi2_recenter->Write();
  // H_Psi3_recenter->Write();
  // H_Psi2_flat->Write();
  // H_Psi3_flat->Write();


  // for(Int_t i=0; i<ngap; i++){
  //   //_west
  //   for(Int_t j=0; j<9; j++){
  //     H_Psi2_west[i][j]->Write();
  //     H_Psi3_west[i][j]->Write();
  //     H_Psi2_recenter_west[i][j]->Write();
  //     H_Psi3_recenter_west[i][j]->Write();
  //     H_Psi2_flat_west[i][j]->Write();
  //     H_Psi3_flat_west[i][j]->Write();
  //     //_east
  //     H_Psi2_east[i][j]->Write();
  //     H_Psi3_east[i][j]->Write();
  //     H_Psi2_recenter_east[i][j]->Write();
  //     H_Psi3_recenter_east[i][j]->Write();
  //     H_Psi2_flat_east[i][j]->Write();
  //     H_Psi3_flat_east[i][j]->Write();
  //   }

  //   P3_Qx2_east[i]->Write();
  //   P3_Qy2_east[i]->Write();
  //   P3_Qx3_east[i]->Write();
  //   P3_Qy3_east[i]->Write();

  //   P3_Qx2_west[i]->Write();
  //   P3_Qy2_west[i]->Write();
  //   P3_Qx3_west[i]->Write();
  //   P3_Qy3_west[i]->Write();

  //   P_square_resolution2[i]->Write();
  //   P_square_resolution3[i]->Write();
  // }
  
  savefile->Write();
  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}
