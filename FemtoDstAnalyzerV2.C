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
void FemtoDstAnalyzerV2(const Char_t *inFile = "AuAu27GeV/AuAu27_ar.list", const Char_t *outFileName = "outFile27.root",const Char_t *flatFileName="Result27Recenter.root", const Char_t *resolutionFileName="Result27Flat.root") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  //Saving our files
  TFile *savefile = new TFile(outFileName, "RECREATE");

  Double_t omega=0;//, Qx2=0, Qy2=0, Qx2_recenter=0, Qy2_recenter=0, Psi2=0, Psi2_recenter=0, delta_Psi2=0, Psi2_flat=0;
  //Double_t Qx3=0, Qy3=0, Qx3_recenter=0, Qy3_recenter=0, Psi3=0, Psi3_recenter=0, delta_Psi3=0, Psi3_flat=0;
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
  Double_t v2=0.0;
  Double_t v3=0.0;

  Double_t Psi2_flat_west[ngap]={0.0}, Psi3_flat_west[ngap]={0.0}, delta_Psi2_west[ngap]={0.0}, delta_Psi3_west[ngap]={0.0};
  Double_t Psi2_flat_east[ngap]={0.0}, Psi3_flat_east[ngap]={0.0}, delta_Psi2_east[ngap]={0.0}, delta_Psi3_east[ngap]={0.0};
  Double_t N_east[ngap]={ 0.0 },   N_west[ngap]={ 0.0 };
  Double_t N=0;
  Int_t q=0;
  Double_t n[ngap]={0.0, 0.075, 0.05, 0.5};  
  Int_t cent[10]={80,70,60,50,40,30,20,10,5,0}; 
  const Char_t *hadrons[]={"Pion","Kaon","Proton"};
  Int_t PID=10;
  Int_t cut_number = 10;

  // Double_t cut_m2_pt[4]={0.7, 0.95, 1.2, 1.5};
  // Double_t cut_m2_mean[3]={0.0194, 0.2443, 0.8816};
  // Double_t cut_m2_sigma_pion[4]={0.005354, 0.01937, 0.03237, 0.04696};
  // Double_t cut_m2_sigma_kaon[4]={0.01218, 0.02349, 0.03465, 0.0564}; 
  // Double_t cut_m2_sigma_proton[4]={0.03302, 0.0398, 0.0528, 0.07038}; 
  // Double_t cut_pt_v_vs_centr_up[3]={1.5, 1.5, 2.5};
  // Double_t cut_pt_v_vs_centr_down[3]={0.2, 0.2, 0.4};
  Double_t cut_m2_pt[4]={0.7, 0.95, 1.2, 1.5};
  Double_t cut_m2_mean[3]={0.0145, 0.2433, 0.8786}; // проверить среднее пионов, оно сдвигается в меньшую сторону при возрастании pt
  Double_t cut_m2_sigma_pion[4]={0.004094, 0.01423, 0.02412, 0.03767};
  Double_t cut_m2_sigma_kaon[4]={0.009838, 0.01741, 0.02779, 0.04554}; 
  Double_t cut_m2_sigma_proton[4]={0.02943, 0.03248, 0.04147, 0.05512}; 
  Double_t cut_pt_v_vs_centr_up[3]={1.5, 1.5, 2.5};
  Double_t cut_pt_v_vs_centr_down[3]={0.2, 0.2, 0.4};

  TProfile2D *P_v2_pt[ngap];
  TProfile2D *P_v3_pt[ngap];
  TProfile2D *P_v2_pt_east[ngap];
  TProfile2D *P_v3_pt_east[ngap];
  TProfile2D *P_v2_pt_west[ngap];
  TProfile2D *P_v3_pt_west[ngap];
  TProfile *P_v2_cent[ngap];
  TProfile *P_v3_cent[ngap];
  TProfile *P_v2_cent_east[ngap];
  TProfile *P_v3_cent_east[ngap];
  TProfile *P_v2_cent_west[ngap];
  TProfile *P_v3_cent_west[ngap];
  
  for(Int_t i=0; i<ngap; i++){
    //Profiles for v2 and v3
    P_v2_pt[i] = new TProfile2D(Form("P_v2_pt_n=%i",i+1),Form("Profile of p_{t} versus v_{2} #Delta#eta-gap=%.2f; p_{t},[Gev/c]; Centrality ", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v3_pt[i] = new TProfile2D(Form("P_v3_pt_n=%i",i+1),Form("Profile of p_{t} versus v_{3} #Delta#eta-gap=%.2f;  p_{t},[Gev/c]; Centrality", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v2_cent[i] = new TProfile(Form("P_v2_cent_n=%i",i+1),Form("Profile of Centrality versus v_{2},#Delta#eta-gap=%.2f; Centrality; v_{2}",2.0*n[i]), 9, -0.5, 8.5);
    P_v3_cent[i] = new TProfile(Form("P_v3_cent_n=%i",i+1),Form("Profile of Centrality versus v_{3},#Delta#eta-gap=%.2f; Centrality; v_{3}",2.0*n[i]), 9, -0.5, 8.5);

    //east
    P_v2_pt_east[i] = new TProfile2D(Form("P_v2_pt_east_n=%i",i+1),Form("Profile of p_{t} versus v_{2} #Delta#eta-gap=%.2f;  p_{t},[Gev/c]; Centrality", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v3_pt_east[i] = new TProfile2D(Form("P_v3_pt_east_n=%i",i+1),Form("Profile of p_{t} versus v_{3} #Delta#eta-gap=%.2f;  p_{t},[Gev/c]; Centrality", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v2_cent_east[i] = new TProfile(Form("P_v2_cent_east_n=%i",i+1),Form("Profile of Centrality versus v_{2},#Delta#eta-gap=%.2f; Centrality; v_{2}",2.0*n[i]), 9, -0.5, 8.5);
    P_v3_cent_east[i] = new TProfile(Form("P_v3_cent_east_n=%i",i+1),Form("Profile of Centrality versus v_{3},#Delta#eta-gap=%.2f; Centrality; v_{3}",2.0*n[i]), 9, -0.5, 8.5);
    //west
    P_v2_pt_west[i] = new TProfile2D(Form("P_v2_pt_west_n=%i",i+1),Form("Profile of p_{t} versus v_{2} #Delta#eta-gap=%.2f;  p_{t},[Gev/c]; Centrality", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v3_pt_west[i] = new TProfile2D(Form("P_v3_pt_west_n=%i",i+1),Form("Profile of p_{t} versus v_{3} #Delta#eta-gap=%.2f;  p_{t},[Gev/c]; Centrality", 2.0*n[i]), 50, 0.2, 5.2, 9, -0.5, 8.5);
    P_v2_cent_west[i] = new TProfile(Form("P_v2_cent_west_n=%i",i+1),Form("Profile of Centrality versus v_{2},#Delta#eta-gap=%.2f; Centrality; v_{2}",2.0*n[i]), 9, -0.5, 8.5);
    P_v3_cent_west[i] = new TProfile(Form("P_v3_cent_west_n=%i",i+1),Form("Profile of Centrality versus v_{3},#Delta#eta-gap=%.2f; Centrality; v_{3}",2.0*n[i]), 9, -0.5, 8.5);
  }



  





  TProfile *p_v2_cent[ngap][3];
  TProfile *p_v2_cent_east[ngap][3];
  TProfile *p_v2_cent_west[ngap][3];
  TProfile *p_v3_cent[ngap][3];
  TProfile *p_v3_cent_east[ngap][3];
  TProfile *p_v3_cent_west[ngap][3];

  TProfile *p_v2_cent_pos[ngap][3];
  TProfile *p_v2_cent_east_pos[ngap][3];
  TProfile *p_v2_cent_west_pos[ngap][3];
  TProfile *p_v3_cent_pos[ngap][3];
  TProfile *p_v3_cent_east_pos[ngap][3];
  TProfile *p_v3_cent_west_pos[ngap][3];

  TProfile *p_v2_cent_neg[ngap][3];
  TProfile *p_v2_cent_east_neg[ngap][3];
  TProfile *p_v2_cent_west_neg[ngap][3];
  TProfile *p_v3_cent_neg[ngap][3];
  TProfile *p_v3_cent_east_neg[ngap][3];
  TProfile *p_v3_cent_west_neg[ngap][3];

  TProfile2D *p_v2_pt[ngap][3];
  TProfile2D *p_v2_pt_east[ngap][3];
  TProfile2D *p_v2_pt_west[ngap][3];
  TProfile2D *p_v3_pt[ngap][3];
  TProfile2D *p_v3_pt_east[ngap][3];
  TProfile2D *p_v3_pt_west[ngap][3];

  TProfile2D *p_v2_pt_pos[ngap][3];
  TProfile2D *p_v2_pt_east_pos[ngap][3];
  TProfile2D *p_v2_pt_west_pos[ngap][3];
  TProfile2D *p_v3_pt_pos[ngap][3];
  TProfile2D *p_v3_pt_east_pos[ngap][3];
  TProfile2D *p_v3_pt_west_pos[ngap][3];

  TProfile2D *p_v2_pt_neg[ngap][3];
  TProfile2D *p_v2_pt_east_neg[ngap][3];
  TProfile2D *p_v2_pt_west_neg[ngap][3];
  TProfile2D *p_v3_pt_neg[ngap][3];
  TProfile2D *p_v3_pt_east_neg[ngap][3];
  TProfile2D *p_v3_pt_west_neg[ngap][3];

  for(Int_t j=0; j<3; j++){
    for(Int_t i=0; i<ngap; i++){

      p_v2_cent[i][j] = new TProfile(Form("v2_cent_prof_n_%i_%s",i+1, hadrons[j]),Form("V_{2}, #Delta#eta-gap=%.2f %s; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_pos[i][j] = new TProfile(Form("v2_cent_prof_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{2}, #Delta#eta-gap=%.2f %s Positive; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_neg[i][j] = new TProfile(Form("v2_cent_prof_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{2}, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v2_cent_east[i][j] = new TProfile(Form("v2_cent_prof_east_n_%i_%s",i+1, hadrons[j]),Form("V_{2} east, #Delta#eta-gap=%.2f %s; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_east_pos[i][j] = new TProfile(Form("v2_cent_prof_east_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{2} east Positive, #Delta#eta-gap=%.2f %s; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_east_neg[i][j] = new TProfile(Form("v2_cent_prof_east_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{2} east, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v2_cent_west[i][j] = new TProfile(Form("v2_cent_prof_west_n_%i_%s",i+1, hadrons[j]),Form("V_{2} west, #Delta#eta-gap=%.2f %s; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_west_pos[i][j] = new TProfile(Form("v2_cent_prof_west_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{2} west Positive, #Delta#eta-gap=%.2f %s; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v2_cent_west_neg[i][j] = new TProfile(Form("v2_cent_prof_west_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{2} west, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{2}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v2_pt[i][j] = new TProfile2D (Form("v2_pt_prof_n_%i_%s",i+1, hadrons[j]),Form("v_{2} #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,50,0.2,5.2);
      p_v2_pt_pos[i][j] = new TProfile2D (Form("v2_pt_prof_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{2} #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v2_pt_neg[i][j] = new TProfile2D (Form("v2_pt_prof_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{2} #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

      p_v2_pt_east[i][j] = new TProfile2D (Form("v2_pt_prof_east_n_%i_%s",i+1, hadrons[j]),Form("v_{2} east #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v2_pt_east_pos[i][j] = new TProfile2D (Form("v2_pt_prof_east_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{2} east #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v2_pt_east_neg[i][j] = new TProfile2D (Form("v2_pt_prof_east_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{2} east #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

      p_v2_pt_west[i][j] = new TProfile2D (Form("v2_pt_prof_west_n_%i_%s",i+1, hadrons[j]),Form("v_{2} west #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v2_pt_west_pos[i][j] = new TProfile2D (Form("v2_pt_prof_west_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{2} west #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v2_pt_west_neg[i][j] = new TProfile2D (Form("v2_pt_prof_west_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{2} west #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

      p_v3_cent[i][j] = new TProfile(Form("v3_cent_prof_n_%i_%s",i+1, hadrons[j]),Form("V_{3}, #Delta#eta-gap=%.2f %s; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_pos[i][j] = new TProfile(Form("v3_cent_prof_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{3}, #Delta#eta-gap=%.2f %s Positive; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_neg[i][j] = new TProfile(Form("v3_cent_prof_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{3}, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v3_cent_east[i][j] = new TProfile(Form("v3_cent_prof_east_n_%i_%s",i+1, hadrons[j]),Form("V_{3} east, #Delta#eta-gap=%.2f %s; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_east_pos[i][j] = new TProfile(Form("v3_cent_prof_east_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{3} east Positive, #Delta#eta-gap=%.2f %s; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_east_neg[i][j] = new TProfile(Form("v3_cent_prof_east_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{3} east, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v3_cent_west[i][j] = new TProfile(Form("v3_cent_prof_west_n_%i_%s",i+1, hadrons[j]),Form("V_{3} west, #Delta#eta-gap=%.2f %s; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_west_pos[i][j] = new TProfile(Form("v3_cent_prof_west_n_%i_%s_pos",i+1, hadrons[j]),Form("V_{3} west Positive, #Delta#eta-gap=%.2f %s; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);
      p_v3_cent_west_neg[i][j] = new TProfile(Form("v3_cent_prof_west_n_%i_%s_neg",i+1, hadrons[j]),Form("V_{3} west, #Delta#eta-gap=%.2f %s Negative; Centrality; V_{3}",2.0*n[i], hadrons[j] ),9,-0.5,8.5);

      p_v3_pt[i][j] = new TProfile2D (Form("v3_pt_prof_n_%i_%s",i+1, hadrons[j]),Form("v_{3}, #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v3_pt_pos[i][j] = new TProfile2D (Form("v3_pt_prof_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{3}, #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v3_pt_neg[i][j] = new TProfile2D (Form("v3_pt_prof_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{3}, #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

      p_v3_pt_east[i][j] = new TProfile2D (Form("v3_pt_prof_east_n_%i_%s",i+1, hadrons[j]),Form("v_{3} east #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v3_pt_east_pos[i][j] = new TProfile2D (Form("v3_pt_prof_east_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{3} east #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);
      p_v3_pt_east_neg[i][j] = new TProfile2D (Form("v3_pt_prof_east_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{3} east #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

      p_v3_pt_west[i][j] = new TProfile2D (Form("v3_pt_prof_west_n_%i_%s",i+1, hadrons[j]),Form("v_{3} west #Delta#eta-gap=%.2f %s; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);  
      p_v3_pt_west_pos[i][j] = new TProfile2D (Form("v3_pt_prof_west_n_%i_%s_pos",i+1, hadrons[j]),Form("v_{3} west #Delta#eta-gap=%.2f %s Positive; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);  
      p_v3_pt_west_neg[i][j] = new TProfile2D (Form("v3_pt_prof_west_n_%i_%s_neg",i+1, hadrons[j]),Form("v_{3} west #Delta#eta-gap=%.2f %s Negative; Centrality; P_{T}[Gev/c]", 2.0*n[i], hadrons[j] ),9,-0.5,8.5,1000,0.2,5.2);

    }
  }
  
  
  
  










 //----------READING FLATENING-------------------------
  TFile *readf = new TFile(flatFileName,"READ");
  
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





  //----------READING RESOLUTION-------------------------  
  TFile *readf2 = new TFile(resolutionFileName,"READ");
  TProfile *P2_square_resolution2[ngap];
  TProfile *P2_square_resolution3[ngap];
  for(Int_t i=0; i<ngap; i++){
    P2_square_resolution2[i] = (TProfile*)readf2 -> Get(Form("P_square_resolution2_n=%i",i+1));
    P2_square_resolution3[i] = (TProfile*)readf2 -> Get(Form("P_square_resolution3_n=%i",i+1));
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
    // Qx2=0; Qy2=0; Psi2=0; Qx3=0; Qy3=0; Psi3=0;
    // Qx2_recenter=0; Qy2_recenter=0; Psi2_recenter=0;
    // Qx3_recenter=0; Qy3_recenter=0; Psi3_recenter=0;
    // delta_Psi3=0; Psi3_flat=0;
    // delta_Psi2=0; Psi2_flat=0;
    v2=0.0;
    v3=0.0;
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
    //N=0;

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

      // Check if track has TOF signal
      if ( femtoTrack->isTofTrack() ) {

      } //if( isTofTrack() )

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    





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

    }


    
    //-------------------Second Track loop for v2 and v3-----------------
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      StFemtoTrack *femtoTrack2 = dst->track(iTrk);
      if (!femtoTrack2) continue;
      if ( !femtoTrack2->isPrimary() ) continue;
      if( (femtoTrack2->dEdx()) == 0 ) continue;
      if( femtoTrack2->gMom().Mag() < 0.1 || femtoTrack2->gDCA(pVtx).Mag() > 2. ) {
         continue;
      }
      //Для индитификации частиц
      if( TMath::Abs( femtoTrack2 -> eta() ) > 1.0 ||
          femtoTrack2 -> nHits() < 15 ||
          femtoTrack2 -> pt() < 0.2 || 
          femtoTrack2 -> p() < 0.15 || 
          femtoTrack2 -> p() > 10.0 ) {
         continue;
      }

      PID=10;
      cut_number=10;

      // // Check if track has TOF signal                   //--------------------вЫключили TPC only---------------------------------------
      // if ( !femtoTrack2->isTofTrack() ){

      //   if( femtoTrack2 -> pt() < 0.55 && TMath::Abs( femtoTrack2->nSigmaPion() ) < 2  ){
      //     PID = 0;
      //   }
      //   if( femtoTrack2 -> pt() < 0.4 && TMath::Abs( femtoTrack2->nSigmaKaon() ) < 2  ){
      //     PID = 1;
      //   }
      //   if( femtoTrack2 -> pt() > 0.2 && femtoTrack2 -> pt() < 0.7 && TMath::Abs( femtoTrack2->nSigmaProton() ) < 2 ){
      //     PID = 2;
      //   }

      // }// if( !femtoTrack2->isTofTrack() )

      // Check if track has TOF signal
      if ( femtoTrack2->isTofTrack() ) {

        if(femtoTrack2->pt() <0.7){
          cut_number = 0;
        }
        if(femtoTrack2->pt() >= 0.7 && femtoTrack2->pt() <0.95){
          cut_number = 1;
        }
        if(femtoTrack2->pt() >= 0.95 && femtoTrack2->pt() <1.2){
          cut_number = 2;
        }
        if(femtoTrack2->pt() >= 1.2 && femtoTrack2->pt() <1.5){
          cut_number = 3;
        }
        if(femtoTrack2->pt() >= 1.5 && femtoTrack2->pt() <3.5){
          cut_number = 4;
        }


        if(cut_number < 4){

          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[0] ) < 2.0*cut_m2_sigma_pion[cut_number] && TMath::Abs( femtoTrack2->nSigmaPion() ) < 3  ){
            PID = 0;
          }
          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[1] ) < 2.0*cut_m2_sigma_kaon[cut_number] && TMath::Abs( femtoTrack2->nSigmaKaon() ) < 3  ){
            PID = 1;
          }
          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[2] ) < 2.0*cut_m2_sigma_proton[cut_number] && TMath::Abs( femtoTrack2->nSigmaProton() ) < 3 ){
            PID = 2;
          }

        }

        if(cut_number == 4){

          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[0] ) < 2.0*cut_m2_sigma_pion[1] && TMath::Abs( femtoTrack2->nSigmaPion() ) < 3  ){
            PID = 0;
          }
          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[1] ) < 2.0*cut_m2_sigma_kaon[1] && TMath::Abs( femtoTrack2->nSigmaKaon() ) < 3  ){
            PID = 1;
          }
          if( TMath::Abs( femtoTrack2->massSqr() - cut_m2_mean[2] ) < 2.0*cut_m2_sigma_proton[1] && TMath::Abs( femtoTrack2->nSigmaProton() ) < 3  ){
            PID = 2;
          }

        }

      } //if( isTofTrack() )





      if ( femtoTrack2->pt()<2){
        omega = femtoTrack2->pt();
      }
      else {
        omega = 2.0;
      }


      if(PID < 5){
        for(Int_t i=0; i<ngap; i++){
          //East TPC
          if(TMath::Abs(femtoTrack2->eta())>n[i] && femtoTrack2->eta() > 0.0 ){
            
            v2 = TMath :: Cos(2.0*(femtoTrack2->phi()-Psi2_west[i])) / TMath :: Sqrt (P2_square_resolution2[i] -> GetBinContent(P2_square_resolution2[i] -> FindBin(event->cent9())));
            v3 = TMath :: Cos(3.0*(femtoTrack2->phi()-Psi3_west[i])) / TMath :: Sqrt (P2_square_resolution3[i] -> GetBinContent(P2_square_resolution3[i] -> FindBin(event->cent9())));
            
            if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){       /// от центральности при 0.2<pt<2.0
              p_v2_cent[i][PID] -> Fill(event->cent9(), v2);
              p_v2_cent_east[i][PID] -> Fill(event->cent9(), v2);
              p_v3_cent[i][PID] -> Fill(event->cent9(), v3);
              p_v3_cent_east[i][PID] -> Fill(event->cent9(), v3);
            }

            p_v2_pt[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
            p_v2_pt_east[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
            p_v3_pt[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
            p_v3_pt_east[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
          
            //------------------------ negative hadrons
            if( femtoTrack2->charge() < 0.0 ){

              if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){       /// от центральности при 0.2<pt<2.0
                p_v2_cent_neg[i][PID] -> Fill(event->cent9(), v2);
                p_v2_cent_east_neg[i][PID] -> Fill(event->cent9(), v2);
                p_v3_cent_neg[i][PID] -> Fill(event->cent9(), v3);
                p_v3_cent_east_neg[i][PID] -> Fill(event->cent9(), v3);
              }

              p_v2_pt_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v2_pt_east_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v3_pt_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
              p_v3_pt_east_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
            }

            //------------------------ \positive hadrons
            if( femtoTrack2->charge() > 0.0 ){
              
              if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){       /// от центральности при 0.2<pt<2.0
                p_v2_cent_pos[i][PID] -> Fill(event->cent9(), v2);
                p_v2_cent_east_pos[i][PID] -> Fill(event->cent9(), v2);
                p_v3_cent_pos[i][PID] -> Fill(event->cent9(), v3);
                p_v3_cent_east_pos[i][PID] -> Fill(event->cent9(), v3);
              }

              p_v2_pt_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v2_pt_east_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v3_pt_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
              p_v3_pt_east_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
            }

          }//east tpc

          //West TPC
          if(TMath::Abs(femtoTrack2->eta())>n[i] && femtoTrack2->eta() < 0.0 ){
            
            v2 = TMath :: Cos(2.0*(femtoTrack2->phi()-Psi2_east[i])) / TMath :: Sqrt (P2_square_resolution2[i] -> GetBinContent(P2_square_resolution2[i] -> FindBin(event->cent9())));
            v3 = TMath :: Cos(3.0*(femtoTrack2->phi()-Psi3_east[i])) / TMath :: Sqrt (P2_square_resolution3[i] -> GetBinContent(P2_square_resolution3[i] -> FindBin(event->cent9())));
            
            if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){     // от центральности при 0.2<pt<2.0
              p_v2_cent[i][PID] -> Fill(event->cent9(), v2);
              p_v2_cent_west[i][PID] -> Fill(event->cent9(), v2);
              p_v3_cent[i][PID] -> Fill(event->cent9(), v3);
              p_v3_cent_west[i][PID] -> Fill(event->cent9(), v3);
            }

            p_v2_pt[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
            p_v2_pt_west[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);          
            p_v3_pt[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
            p_v3_pt_west[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
            
            // -----------------------------------negative hadrons
            if( femtoTrack2->charge() < 0 ){

              if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){     // от центральности при 0.2<pt<2.0
                p_v2_cent_neg[i][PID] -> Fill(event->cent9(), v2);
                p_v2_cent_west_neg[i][PID] -> Fill(event->cent9(), v2);
                p_v3_cent_neg[i][PID] -> Fill(event->cent9(), v3);
                p_v3_cent_west_neg[i][PID] -> Fill(event->cent9(), v3);
              }

              p_v2_pt_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v2_pt_west_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);          
              p_v3_pt_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
              p_v3_pt_west_neg[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);

            }

            //------------------------------------positive hadrons
            if( femtoTrack2->charge() > 0 ){

              if(femtoTrack2->pt() < cut_pt_v_vs_centr_up[PID] && femtoTrack2->pt() > cut_pt_v_vs_centr_down[PID] ){     // от центральности при 0.2<pt<2.0
                p_v2_cent_pos[i][PID] -> Fill(event->cent9(), v2);
                p_v2_cent_west_pos[i][PID] -> Fill(event->cent9(), v2);
                p_v3_cent_pos[i][PID] -> Fill(event->cent9(), v3);
                p_v3_cent_west_pos[i][PID] -> Fill(event->cent9(), v3);
              }

              p_v2_pt_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);
              p_v2_pt_west_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v2, omega);          
              p_v3_pt_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);
              p_v3_pt_west_pos[i][PID] -> Fill(event->cent9(),femtoTrack2->pt(), v3, omega);

            }
          
          }//west tpc
        }//gap

      }// if( PID < 5)
    }//end second track loop

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


  
  savefile->Write();
  savefile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
}