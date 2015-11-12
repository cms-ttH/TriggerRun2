#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;

double reweightPU( int nPU );
float DeltaR(float eta1,float phi1,float eta2,float phi2);

// ------------ csv applying functions -------------
void fillCSVhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

int PtBinsHF_ = 5;

// CSV reweighting
TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];

//*****************************************************************************

void ttHbb_data2mc_treeReader( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, int ttCat_=-1 ) {

  std::string inputFileHF = "data/csv_rwt_hf_IT_FlatSF_2015_07_27.root";
  std::string inputFileLF = "data/csv_rwt_lf_IT_FlatSF_2015_07_27.root";

  // TFile* f_CSVwgt_HF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/TriggerRun2/TriggerAnalyzer/" + inputFileHF).c_str());
  // TFile* f_CSVwgt_LF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/TriggerRun2/TriggerAnalyzer/" + inputFileLF).c_str());

  // TFile* f_CSVwgt_HF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_5/src/ttH-LeptonPlusJets/AnalysisCode/csv_rwt_fit_hf_v2_final_2015_11_03.root");
  // TFile* f_CSVwgt_LF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_5/src/ttH-LeptonPlusJets/AnalysisCode/csv_rwt_fit_lf_v2_final_2015_11_03.root");

  TFile* f_CSVwgt_HF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_12/src/csvReweightingRun2/csvTreeMaker/csv_rwt_fit_hf_v3_final_2015_11_12.root");
  TFile* f_CSVwgt_LF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_12/src/csvReweightingRun2/csvTreeMaker/csv_rwt_fit_lf_v3_final_2015_11_12.root");


  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

  ////
  std::cout << "   ===> load the root files! " << std::endl;

  std::string sampleType = ( insample>=0 ) ? "mc" : "data";
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
  double mySample_xSec_ = 1.;
  double mySample_nGen_ = 1.;
  std::string mySample_sampleName_ = "delete";
  std::string mySample_inputDir_ = "";
  std::vector<std::string> mySample_inputDirs_;
  if( insample==2500 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 115091972;//25357774;//25446993;
    mySample_sampleName_ = "ttbar";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_triggerTree_v1/151015_183717/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1_triggerTree_v1/151015_183747/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb";
    }
  }
  else if( insample==2300 ){
    mySample_xSec_ = 6025.2; 
    mySample_nGen_ = 19310834;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/151015_183818/0000/");
  }
  else if( insample==2310 ){
    mySample_xSec_ = 22635.09; 
    mySample_nGen_ = 22217467;//19925500;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M10to50";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_183844/0000/");
  }
  else if( insample==2400 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 16518218;
    mySample_sampleName_ = "WJets";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184003/0000/");
  }
  else if( insample==2524 ){
    mySample_xSec_ = 0.435;  
    mySample_nGen_ = 430330;
    mySample_sampleName_ = "ttW_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184416/0000/");
  }
  else if( insample==2525 ){
    mySample_xSec_ = 0.21;  
    mySample_nGen_ = 129850;
    mySample_sampleName_ = "ttW_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184445/0000/");
  }
  else if( insample==2523 ){
    mySample_xSec_ = 0.611;  
    mySample_nGen_ = 351398;
    mySample_sampleName_ = "ttZ_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184313/0000/");
  }
  else if( insample==2522 ){
    mySample_xSec_ = 0.263;  
    mySample_nGen_ = 184990;
    mySample_sampleName_ = "ttZ_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184349/0000/");
  }
  else if( insample==2510 ){
    mySample_xSec_ = 136.02;  
    mySample_nGen_ = 3299800;
    mySample_sampleName_ = "st_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184110/0000/");
  }
  else if( insample==2511 ){
    mySample_xSec_ = 80.95;  
    mySample_nGen_ = 1695400;
    mySample_sampleName_ = "stbar_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184034/0000/");
  }
  else if( insample==2512 ){
    mySample_xSec_ = 35.9;  
    mySample_nGen_ = 995600;
    mySample_sampleName_ = "st_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184134/0000/");
  }
  else if( insample==2513 ){
    mySample_xSec_ = 35.9;  
    mySample_nGen_ = 1000000;
    mySample_sampleName_ = "stbar_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184202/0000/");
  }
  else if( insample==2514 ){
    mySample_xSec_ = 10.32;  
    mySample_nGen_ = 613384;
    mySample_sampleName_ = "st_schan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184243/0000/");
  }
  else if( insample==2424 ){
    mySample_xSec_ = 118.7;  
    mySample_nGen_ = 994416;
    mySample_sampleName_ = "WW";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184512/0000/");
  }
  else if( insample==2423 ){
    mySample_xSec_ = 44.9;  
    mySample_nGen_ = 991232;
    mySample_sampleName_ = "WZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184539/0000/");
  }
  else if( insample==2323 ){
    mySample_xSec_ = 15.4;  
    mySample_nGen_ = 996168;
    mySample_sampleName_ = "ZZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/151015_184614/0000/");
  }
  else if( insample==9125 ){
    mySample_xSec_ = 0.2934;// YR3 * BR(all)  
    mySample_nGen_ = 3933404;//199000;
    mySample_sampleName_ = "ttHTobb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_183910/0000/");
  }
  else if( insample==8125 ){
    mySample_xSec_ = 0.2151;// YR3 * BR(all)  
    mySample_nGen_ = 3796398;//199000;
    mySample_sampleName_ = "ttHnonbb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_triggerTree_v1/151015_183937/0000/");
  }
  else if( insample==-13 ){
    mySample_sampleName_ = "SingleElectron_Run2015D_254231_258158";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v3_GoldenJSON_254231_258158_triggerTree_v1/151012_183037/0000/");
  }
  else if( insample==-11 ){
    mySample_sampleName_ = "SingleMuon_Run2015D_254231_258158";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleMuon/Run2015D-PromptReco-v3_GoldenJSON_254231_258158_triggerTree_v1/151015_184732/0000/");
  }




  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/ttHbb_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());

  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    //chain->Add(treefilename.c_str());
  }
  chain->Add("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1_triggerTree_v1/151015_183747/0000/trigger_analyzer_223.root");

  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  //////////////////////////////////////////////////////////////////////////
  ///  Histogram making
  //////////////////////////////////////////////////////////////////////////


  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();


  bool verbose = false;

  //////////////////////////////////////////////////////////////////////////
  ///  Histograms
  //////////////////////////////////////////////////////////////////////////

  double window = 10.;
  double pdgZmass = 91.1876;
  double MinMass  = pdgZmass - window;
  double MaxMass  = pdgZmass + window;

  TH1::SetDefaultSumw2();

  TH1D* h_numEvents = new TH1D("h_numEvents",";Number of events", 2, 0, 2 );
  TH1D* h_numEvents_wgt = new TH1D("h_numEvents_wgt",";Predicted number of events", 7, 0, 7 );

  TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );
  TH1D* h_numPVs_PUwgt = new TH1D("h_numPVs_PUwgt",";Number of PVs", 50, 0-0.5, 50-0.5 );


  TH1D* h_additionalJetEventId = new TH1D("h_additionalJetEventId",";additionalJetEventId", 201, -100-0.5, 101-0.5 );



  int NumCuts = 5;
  TH1D* h_mu_event_selection  = new TH1D("h_mu_event_selection",";cut", NumCuts, 0, NumCuts );
  h_mu_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_mu_event_selection->GetXaxis()->SetBinLabel(2,"HLT Mu");
  h_mu_event_selection->GetXaxis()->SetBinLabel(3,"==1 muon");
  h_mu_event_selection->GetXaxis()->SetBinLabel(4,">=4 jets");
  h_mu_event_selection->GetXaxis()->SetBinLabel(5,">=2 b-jets");

  TH1D* h_ele_event_selection  = new TH1D("h_ele_event_selection",";cut", NumCuts, 0, NumCuts );
  h_ele_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_ele_event_selection->GetXaxis()->SetBinLabel(2,"HLT Ele");
  h_ele_event_selection->GetXaxis()->SetBinLabel(3,"==1 electron");
  h_ele_event_selection->GetXaxis()->SetBinLabel(4,">=4 jets");
  h_ele_event_selection->GetXaxis()->SetBinLabel(5,">=2 b-jets");

  TH1D* h_event_selection  = new TH1D("h_event_selection",";cut", NumCuts, 0, NumCuts );
  h_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_event_selection->GetXaxis()->SetBinLabel(2,"HLT");
  h_event_selection->GetXaxis()->SetBinLabel(3,"==1 lepton");
  h_event_selection->GetXaxis()->SetBinLabel(4,">=4 jets");
  h_event_selection->GetXaxis()->SetBinLabel(5,">=2 b-jets");


  TH1D* h_deltaR_jet_lep = new TH1D("h_deltaR_jet_lep",";#DeltaR(jet,lep)", 61, 0., 6.1 );

  TH1D* h_numJet = new TH1D("h_numJet",";Number of Jets", 5, 4-0.5, 9-0.5 );
  TH1D* h_numBtag = new TH1D("h_numBtag",";Number of b-tagged Jets", 4, 2-0.5, 6-0.5 );


  double HTmax = 1000.;
  int numHTbins = 1000;
  double L1HTmax = HTmax;
  int numL1HTTbins = numHTbins;


  std::vector<std::string> cat_labels;
  cat_labels.push_back("incl4j2t");
  cat_labels.push_back("4j2t");
  cat_labels.push_back("5j2t");
  cat_labels.push_back("6j2t");
  cat_labels.push_back("4j3t");
  cat_labels.push_back("5j3t");
  cat_labels.push_back("6j3t");
  cat_labels.push_back("4j4t");
  cat_labels.push_back("5j4t");
  cat_labels.push_back("6j4t");

  int NumCat = int(cat_labels.size());

  TH1D* h_category_yield = new TH1D("h_category_yield", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e = new TH1D("h_category_yield_1e", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1m = new TH1D("h_category_yield_1m", ";category", NumCat, 0, NumCat );

  for( int iCat=0; iCat<NumCat; iCat++ ){
    h_category_yield->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_1e->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_1m->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
  }


  std::vector<int> systematics;
  systematics.push_back(0);
  systematics.push_back(7);
  systematics.push_back(8);
  systematics.push_back(9);
  systematics.push_back(10);
  systematics.push_back(11);
  systematics.push_back(12);
  systematics.push_back(13);
  systematics.push_back(14);
  systematics.push_back(15);
  systematics.push_back(16);
  systematics.push_back(17);
  systematics.push_back(18);
  systematics.push_back(19);
  systematics.push_back(20);
  systematics.push_back(21);
  systematics.push_back(22);
  systematics.push_back(23);
  systematics.push_back(24);

  int numSys = int(systematics.size());


  TH1D* h_numEvents_Sys = new TH1D("h_numEvents_Sys",";systematic;predicted number of events", systematics[numSys-1]+1, 0, systematics[numSys-1]+1 );
  TH1D* h_numEvents_4j2t_Sys = new TH1D("h_numEvents_4j2t_Sys",";systematic;predicted number of events", systematics[numSys-1]+1, 0, systematics[numSys-1]+1 );

  TH1D* h_numEvents_Sys_PU = new TH1D("h_numEvents_Sys_PU",";systematic;predicted number of events", systematics[numSys-1]+1, 0, systematics[numSys-1]+1 );

  TH1D* h_wgt_csv_Sys[numSys];
  TH1D* h_wgt_csv_hf_Sys[numSys];
  TH1D* h_wgt_csv_lf_Sys[numSys];
  TH1D* h_wgt_csv_cf_Sys[numSys];

  TH1D* h_wgt_csv_4j2t_Sys[numSys];
  TH1D* h_wgt_csv_hf_4j2t_Sys[numSys];
  TH1D* h_wgt_csv_lf_4j2t_Sys[numSys];
  TH1D* h_wgt_csv_cf_4j2t_Sys[numSys];
  for( int iSys=0; iSys<numSys; iSys++ ){
    int useSys = systematics[iSys];

    h_wgt_csv_Sys[iSys] = new TH1D(Form("h_wgt_csv_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_Sys_%d",useSys),";weight", 300, 0., 3. );

    h_wgt_csv_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
  }


  //
  // Lepton plots
  //
  TH1D* h_lepton_pt[NumCat];
  TH1D* h_lepton_phi[NumCat];
  TH1D* h_lepton_eta[NumCat];

  TH1D* h_electron_pt[NumCat];
  TH1D* h_electron_phi[NumCat];
  TH1D* h_electron_eta[NumCat];

  TH1D* h_muon_pt[NumCat];
  TH1D* h_muon_phi[NumCat];
  TH1D* h_muon_eta[NumCat];


  // 
  // Jet plots
  //
  TH1D* h_jet_pt[NumCat];
  TH1D* h_jet_eta[NumCat];
  TH1D* h_jet_phi[NumCat];
  TH1D* h_jet_csv[NumCat];

  //
  // Energy sum plots
  //
  TH1D* h_met_pt[NumCat];
  TH1D* h_met_phi[NumCat];
  TH1D* h_mht_pt[NumCat];
  TH1D* h_mht_phi[NumCat];

  TH1D* h_HT[NumCat];


  //// Histograming 
  double metmax   = 500.;
  double lepPtMax = 300.;
  double jetptmax = 500.;
  double htmax    = 2000.;
  int NcsvBins = 44;

  int NmetBins   = int( metmax/20. + 0.0001 );
  int NlepPtBins = int( lepPtMax/10. + 0.0001 );
  int NjetptBins = int( jetptmax/10. + 0.0001 );
  int NhtBins    = int( htmax/50. + 0.0001 );

  for( int c=0; c<NumCat; c++ ){ 
    std::string suffix = "_" + cat_labels[c];
    std::string cat_suffix = "_" + cat_labels[c];

    h_lepton_pt[c] = new TH1D((std::string("h_lepton_pt" + suffix)).c_str(),";lepton p_{T}", NlepPtBins, 0, lepPtMax );
    h_lepton_phi[c] = new TH1D((std::string("h_lepton_phi" + suffix)).c_str(),";lepton #phi", 34, -3.4, 3.4 );
    h_lepton_eta[c] = new TH1D((std::string("h_lepton_eta" + suffix)).c_str(),";lepton #eta", 25, -2.5, 2.5 );

    h_electron_pt[c] = new TH1D((std::string("h_electron_pt" + suffix)).c_str(),";electron p_{T}", NlepPtBins, 0, lepPtMax );
    h_electron_phi[c] = new TH1D((std::string("h_electron_phi" + suffix)).c_str(),";electron #phi", 34, -3.4, 3.4 );
    h_electron_eta[c] = new TH1D((std::string("h_electron_eta" + suffix)).c_str(),";electron #eta", 25, -2.5, 2.5 );

    h_muon_pt[c] = new TH1D((std::string("h_muon_pt" + suffix)).c_str(),";muon p_{T}", NlepPtBins, 0, lepPtMax );
    h_muon_phi[c] = new TH1D((std::string("h_muon_phi" + suffix)).c_str(),";muon #phi", 34, -3.4, 3.4 );
    h_muon_eta[c] = new TH1D((std::string("h_muon_eta" + suffix)).c_str(),";muon #eta", 25, -2.5, 2.5 );

    h_jet_pt[c] = new TH1D((std::string("h_jet_pt" + suffix)).c_str(),";jet p_{T}", NjetptBins, 0, jetptmax );
    h_jet_eta[c] = new TH1D((std::string("h_jet_eta" + suffix)).c_str(),";jet #eta", 25, -2.5, 2.5 );
    h_jet_phi[c] = new TH1D((std::string("h_jet_phi" + suffix)).c_str(),";jet #phi", 34, -3.4, 3.4 );
    h_jet_csv[c] = new TH1D((std::string("h_jet_csv" + suffix)).c_str(),";jet CSV", NcsvBins, -1.1, 1.1 );

    h_met_pt[c]  = new TH1D((std::string("h_met_pt" + suffix)).c_str(),";MET p_{T}", NmetBins, 0, metmax );
    h_met_phi[c] = new TH1D((std::string("h_met_phi" + suffix)).c_str(),";MET #phi", 34, -3.4, 3.4 );
    h_mht_pt[c]  = new TH1D((std::string("h_mht_pt" + suffix)).c_str(),";MHT p_{T}", NmetBins, 0, metmax );
    h_mht_phi[c] = new TH1D((std::string("h_mht_phi" + suffix)).c_str(),";MHT #phi", 34, -3.4, 3.4 );

    h_HT[c] = new TH1D((std::string("h_HT" + suffix)).c_str(),";H_{T} (jets)", NhtBins, 0, htmax );
  }


  int NumHists = 6;
  int NumSysC = 5;
  TH1D* h_c_jet_csv[NumHists];
  TH1D* h_c_jet_csv_CharmCSVSF[NumSysC][NumHists];

  int NfullcsvBins = 1000;
  for( int iPt=0; iPt<NumHists; iPt++ ){
    h_c_jet_csv[iPt] = new TH1D(Form("h_c_jet_csv_Pt%d_Eta%d", iPt, 0),";jet CSV", NfullcsvBins, -0.1, 1.1 );
    for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
      h_c_jet_csv_CharmCSVSF[iSysC][iPt] = new TH1D(Form("h_c_jet_csv_CharmCSVSF_Pt%d_Eta%d_SysC%d", iPt, 0, iSysC),";jet CSV", NfullcsvBins, -0.1, 1.1 );
    }
  }


  //////////////////////////////////////////////////////////////////////////
  /////
  //////////////////////////////////////////////////////////////////////////

  std::string ele_path_name = ( insample<0 ) ? "HLT_Ele27_eta2p1_WPLoose_Gsf_v" : "HLT_Ele27_WP85_Gsf_v";
  std::string mu_path_name  = ( insample<0 ) ? "HLT_IsoMu20_v" : "HLT_IsoMu20_v";

  int numEvents_all = 0;
  int numEvents_wgt_gen = 0;
  double numEvents_wgt_lumi = 0;
  double numEvents_wgt = 0;

  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;

  int cnt = 0;
  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t
    cnt++;
    if( ievt<firstEvent ) continue;
    if( ievt==lastEvent ) break;

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%10000==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					     << int(double(ievt-firstEvent)/double(NeventsPerJob)*100) << "% done" << std::endl;

    //if( ievt==(maxNentries+1) ) break;
    if( ievt==(maxNentries+1) && ievt!=0 ) break;

    chain->GetEntry(ievt);


    int additionalJetEventId = -99;
    if( insample==2500 ){
      additionalJetEventId = eve->additionalJetEventId_;

      bool keepTTbarEvent = true;
      if( ttCat_>=0 ){
	if( ttCat_==0 && additionalJetEventId==0 ) keepTTbarEvent = true;
	else if( ttCat_==1 && (additionalJetEventId>=41 && additionalJetEventId<=45) ) keepTTbarEvent = true;
	else if( ttCat_==2 && additionalJetEventId==51 ) keepTTbarEvent = true;
	else if( ttCat_==3 && additionalJetEventId==52 ) keepTTbarEvent = true;
	else if( ttCat_==4 && (additionalJetEventId>=53 && additionalJetEventId<=55) ) keepTTbarEvent = true;
	else keepTTbarEvent = false;
      }

      // if( keepTTbarEvent==false && ttCat_>=0 ){
      // 	std::cout << " ERROR!! keepTTbarEvent == false! additionalJetEventId = " << additionalJetEventId << " and ttCat_ = " << ttCat_ << std::endl;
      // }

      if( !keepTTbarEvent ) continue;
      //std::cout << " I have made it past the gate! " << std::endl;
    }

    h_additionalJetEventId->Fill(additionalJetEventId);

    // if( ievt<100 ){
    //   printf(" Event %3lld: additionalJetEventId = %3d \n", ievt, additionalJetEventId);
    // }

    std::vector<double> rejetPts;
    std::vector<double> rejetEtas;
    std::vector<double> rejetCSVs;
    std::vector<int>    rejetFlavors;

    int temp_njet = 0;
    int temp_nbtag = 0;

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      int flavor = eve->jet_flavor_[iJet];
      double csv =  eve->jet_csv_[iJet];
      if( csv<0 && csv>-9 ) csv = -0.2;
      else if( csv < -5 )   csv = -0.4;

      if( csv > 1.0 ) csv = 1.0;

      if( abs(flavor)==4 && fabs(eta)<2.4 && pt>20. ){
	int iPt = -1;
	if (pt >=19.99 && pt<30) iPt = 0;
	else if (pt >=30 && pt<40) iPt = 1;
	else if (pt >=40 && pt<60) iPt = 2;
	else if (pt >=60 && pt<100) iPt = 3;
	else if (pt >=100)          iPt = 4;
	//else if (pt >=100 && pt<160) iPt = 4;
	//else if (pt >=160 && pt<10000) iPt = 5;

	if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

	h_c_jet_csv[iPt]->Fill(csv);
	for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
	  int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	  double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	  h_c_jet_csv_CharmCSVSF[iSysC][iPt]->Fill(csv,iCSVWgtC);
	}
      }

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      temp_njet++;
      if( csv>0.89 ) temp_nbtag++;

      rejetPts.push_back(pt);
      rejetEtas.push_back(eta);
      rejetCSVs.push_back(csv);
      rejetFlavors.push_back(flavor);
    }


    int numPVs = eve->numPVs_;

    double scalePU = ( insample < 0 ) ? 1. : reweightPU(numPVs);

    double wgt_csv = 1;
    for( int useSys=0; useSys<int(systematics.size()); useSys++ ){
      int mySys = systematics[useSys];

      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
      double temp_wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(rejetPts, rejetEtas, rejetCSVs, rejetFlavors, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      if( useSys==0 ) wgt_csv = temp_wgt_csv;

      h_numEvents_Sys->Fill(0.5+mySys,temp_wgt_csv);
      h_numEvents_Sys_PU->Fill(0.5+mySys,temp_wgt_csv*scalePU);

      h_wgt_csv_Sys[useSys]->Fill(temp_wgt_csv,scalePU);
      h_wgt_csv_hf_Sys[useSys]->Fill(wgt_csv_hf,scalePU);
      h_wgt_csv_lf_Sys[useSys]->Fill(wgt_csv_lf,scalePU);
      h_wgt_csv_cf_Sys[useSys]->Fill(wgt_csv_cf,scalePU);


      if( temp_njet>=4 && temp_nbtag>=2 ){
	h_numEvents_4j2t_Sys->Fill(0.5+mySys,temp_wgt_csv*scalePU);
	h_wgt_csv_4j2t_Sys[useSys]->Fill(temp_wgt_csv,scalePU);
	h_wgt_csv_hf_4j2t_Sys[useSys]->Fill(wgt_csv_hf,scalePU);
	h_wgt_csv_lf_4j2t_Sys[useSys]->Fill(wgt_csv_lf,scalePU);
	h_wgt_csv_cf_4j2t_Sys[useSys]->Fill(wgt_csv_cf,scalePU);
      }
    }

    if( insample<0 ) wgt_csv = 1.0;

    int run = eve->run_;

    //// --------- various weights: PU, topPt, triggerSF, leptonSF...
    // double  wgt_topPtSF = eve->wgt_topPt_;
    double Xsec = mySample_xSec_;//eve->wgt_xs_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    if( maxNentries==1000000 && insample==2300 ) nGen = 670121;
    double lumi = ( intLumi > 0 ) ? intLumi : 10000 ;

    double wgt_gen = ( insample > 0 ) ? eve->wgt_generator_ : 1;
    wgt_gen = ( wgt_gen > 0 ) ? 1. : -1.;
    double wgt_lumi = ( insample > 0 ) ? lumi * (Xsec/nGen) : 1;//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights

    double wgt = wgt_gen * wgt_lumi * wgt_csv;

    h_numEvents->Fill(0.5,1.);
    h_numEvents->Fill(1.5,wgt_gen);

    h_numEvents_wgt->Fill(0.5,1.);
    h_numEvents_wgt->Fill(1.5,wgt_gen);
    h_numEvents_wgt->Fill(2.5,wgt_lumi);
    h_numEvents_wgt->Fill(3.5,wgt_gen * wgt_lumi);
    h_numEvents_wgt->Fill(4.5,wgt);


    numEvents_all += 1.;
    numEvents_wgt_gen += wgt_gen;
    numEvents_wgt_lumi += wgt_lumi;
    numEvents_wgt += wgt;


    ///////////////////
    ////// selections
    ///////////////////


    h_event_selection->Fill(0.5, wgt); // all
    h_mu_event_selection->Fill(0.5, wgt); // all
    h_ele_event_selection->Fill(0.5, wgt); // all

    // bool pass_trigger = false;
    // if( insample<0 ) pass_trigger = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    // else             pass_trigger = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );
    

    bool pass_trigger_ele = false;
    bool pass_trigger_mu  = false;

    vint hlt_accept = eve->hlt_accept_;
    vstring hlt_name = eve->hlt_name_;

    bool found_trigger_ele = false;
    bool found_trigger_mu  = false;

    int Nhlt = int( hlt_accept.size() );
    for( int iHLT=0; iHLT<Nhlt; iHLT++ ){
      if( found_trigger_ele && found_trigger_mu ) break;

      std::string name = hlt_name[iHLT];
      int accept = hlt_accept[iHLT];

      if( name.find(ele_path_name)!=std::string::npos ){
	found_trigger_ele = true;
	if( accept ) pass_trigger_ele = true;
	else         pass_trigger_ele = false;
      }
      else if( name.find(mu_path_name)!=std::string::npos ){
	found_trigger_mu = true;
	if( accept ) pass_trigger_mu = true;
	else         pass_trigger_mu = false;
      }
    }


    if( insample<0 ){
      pass_trigger_ele = ( insample==-13 && pass_trigger_ele );
      pass_trigger_mu  = ( insample==-11 && pass_trigger_mu  );
    }

    if( pass_trigger_mu || pass_trigger_ele ) h_event_selection->Fill(0.5+1, wgt); // all
    if( pass_trigger_mu )  h_mu_event_selection->Fill(0.5+1, wgt); // all
    if( pass_trigger_ele ) h_ele_event_selection->Fill(0.5+1, wgt); // all

    // Pass Trigger selection
    if( !(pass_trigger_mu || pass_trigger_ele) ) continue;

    h_numEvents_wgt->Fill(5.5,wgt);


    h_numPVs->Fill(numPVs,wgt);


    wgt *= scalePU;

    h_numPVs_PUwgt->Fill(numPVs,wgt);

    h_numEvents_wgt->Fill(6.5,wgt);

    //
    // Lepton selection
    //

    std::vector<int> ind_ele;
    std::vector<int> ind_mu;
    std::vector<int> ind_ele_loose;
    std::vector<int> ind_mu_loose;
    for( int iLep=0; iLep<int(eve->lepton_pt_.size()); iLep++ ){
      bool isMuon = ( eve->lepton_isMuon_[iLep]==1 );

      bool isSpring15M = ( eve->lepton_isSpring15M_[iLep]==1 );
      bool isTight = ( eve->lepton_isTight_[iLep]==1 );
      bool isCrack = ( eve->lepton_inCrack_[iLep]==1 );

      double pt  = eve->lepton_pt_[iLep];
      double eta = eve->lepton_eta_[iLep];

      if( isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isTight ) ind_mu.push_back(iLep);
	if( pt > 20 && abs(eta)<2.4 && isTight ) ind_mu_loose.push_back(iLep);
      }

      if( !isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isSpring15M && !isCrack ) ind_ele.push_back(iLep);
	if( pt > 20 && abs(eta)<2.5 && isSpring15M && !isCrack ) ind_ele_loose.push_back(iLep);
      }
    }

    int numEle = int( ind_ele.size() );
    int numMu  = int( ind_mu.size() );

    int numEle_loose = int( ind_ele_loose.size() );
    int numMu_loose  = int( ind_mu_loose.size() );

    bool oneEle = ( numEle==1 && numEle_loose==1 && numMu_loose==0 );
    bool oneMu  = ( numMu==1 && numMu_loose==1 && numEle_loose==0 );

    oneEle = ( oneEle && pass_trigger_ele );
    oneMu  = ( oneMu  && pass_trigger_mu );

    bool oneLep = ( oneEle || oneMu );

    if( oneLep ) h_event_selection->Fill(0.5+2, wgt); // all
    if( oneMu )  h_mu_event_selection->Fill(0.5+2, wgt); // all
    if( oneEle ) h_ele_event_selection->Fill(0.5+2, wgt); // all

    // Require exactly one lepton
    if( !oneLep ) continue;

    int lepInd = ( numEle>0 ) ? ind_ele[0] : ind_mu[0];

    TLorentzVector myLep;
    myLep.SetPtEtaPhiE( eve->lepton_pt_[lepInd], eve->lepton_eta_[lepInd], eve->lepton_phi_[lepInd], eve->lepton_energy_[lepInd] );


    vdouble jetPts;
    int numJet = 0;
    int numBtag = 0;
    double HT30=0;
    TLorentzVector sumJet;
    bool firstJet = false;
    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      TLorentzVector myJet;
      myJet.SetPtEtaPhiE( eve->jet_pt_[iJet], eve->jet_eta_[iJet], eve->jet_phi_[iJet], eve->jet_energy_[iJet] );

      if( !firstJet ) sumJet = myJet;
      else            sumJet += myJet;

      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      double phi = eve->jet_phi_[iJet];

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      double dR = myJet.DeltaR(myLep);
      h_deltaR_jet_lep->Fill(dR,wgt);

      if( dR < 0.4 ) continue;

      double csv =  eve->jet_csv_[iJet];
      if( csv<0 && csv>-9 ) csv = -0.2;
      else if( csv < -5 )   csv = -0.4;

      if( csv > 1.0 ) csv = 1.0;

      HT30 += pt;

      jetPts.push_back(pt);

      numJet++;
      if( csv > 0.890 ) numBtag++;
    }


    if( !(numJet>=4 && numBtag>=2) ) continue;


    // Require at least four jets
    if( !(numJet>=4) ) continue;

    if( oneLep ) h_event_selection->Fill(0.5+3, wgt); // all
    if( oneMu )  h_mu_event_selection->Fill(0.5+3, wgt); // all
    if( oneEle ) h_ele_event_selection->Fill(0.5+3, wgt); // all

    // Require at least two b-tags
    if( !(numBtag>=2) ) continue;

    if( oneLep ) h_event_selection->Fill(0.5+4, wgt); // all
    if( oneMu )  h_mu_event_selection->Fill(0.5+4, wgt); // all
    if( oneEle ) h_ele_event_selection->Fill(0.5+4, wgt); // all



    ////////////



    // std::sort(jetPts.begin(), jetPts.end());
    // std::reverse(jetPts.begin(), jetPts.end());

    // if( jetPts.size()>=1 ) h_jet_1_pt->Fill(jetPts[0],wgt);
    // if( jetPts.size()>=2 ) h_jet_2_pt->Fill(jetPts[1],wgt);
    // if( jetPts.size()>=3 ) h_jet_3_pt->Fill(jetPts[2],wgt);
    // if( jetPts.size()>=4 ) h_jet_4_pt->Fill(jetPts[3],wgt);

    h_numJet->Fill(numJet,wgt);
    h_numBtag->Fill(numBtag,wgt);


    int this_category = -1;
    if( numJet==4 && numBtag==2) this_category=1;
    if( numJet==5 && numBtag==2) this_category=2;
    if( numJet>=6 && numBtag==2) this_category=3;	
    if( numJet==4 && numBtag==3) this_category=4;
    if( numJet==5 && numBtag==3) this_category=5;
    if( numJet>=6 && numBtag==3) this_category=6;
    if( numJet==4 && numBtag>=4) this_category=7;
    if( numJet==5 && numBtag>=4) this_category=8;
    if( numJet>=6 && numBtag>=4) this_category=9;

    h_category_yield->Fill(0.5,wgt);
    h_category_yield->Fill(this_category,wgt);

    if( oneEle ){
      h_category_yield_1e->Fill(0.5,wgt);
      h_category_yield_1e->Fill(this_category,wgt);
    }
    else if( oneMu ){
      h_category_yield_1m->Fill(0.5,wgt);
      h_category_yield_1m->Fill(this_category,wgt);
    }


    double lepPt = myLep.Pt();
    double lepEta = myLep.Eta();
    double lepPhi = myLep.Phi();

    h_lepton_pt[0]->Fill(lepPt,wgt);
    h_lepton_eta[0]->Fill(lepEta,wgt);
    h_lepton_phi[0]->Fill(lepPhi,wgt);

    h_lepton_pt[this_category]->Fill(lepPt,wgt);
    h_lepton_eta[this_category]->Fill(lepEta,wgt);
    h_lepton_phi[this_category]->Fill(lepPhi,wgt);

    if( oneEle ){
      h_electron_pt[0]->Fill(lepPt,wgt);
      h_electron_eta[0]->Fill(lepEta,wgt);
      h_electron_phi[0]->Fill(lepPhi,wgt);

      h_electron_pt[this_category]->Fill(lepPt,wgt);
      h_electron_eta[this_category]->Fill(lepEta,wgt);
      h_electron_phi[this_category]->Fill(lepPhi,wgt);
    }
    else if( oneMu ){
      h_muon_pt[0]->Fill(lepPt,wgt);
      h_muon_eta[0]->Fill(lepEta,wgt);
      h_muon_phi[0]->Fill(lepPhi,wgt);

      h_muon_pt[this_category]->Fill(lepPt,wgt);
      h_muon_eta[this_category]->Fill(lepEta,wgt);
      h_muon_phi[this_category]->Fill(lepPhi,wgt);
    }



    double met = eve->pfMET_pt_;
    double met_phi = eve->pfMET_phi_;
    double mht = sumJet.Pt();
    double mht_phi = sumJet.Phi();


    h_met_pt[0]->Fill(met,wgt);
    h_met_pt[this_category]->Fill(met,wgt);

    h_met_phi[0]->Fill(met_phi,wgt);
    h_met_phi[this_category]->Fill(met_phi,wgt);

    h_mht_pt[0]->Fill(mht,wgt);
    h_mht_pt[this_category]->Fill(mht,wgt);

    h_mht_phi[0]->Fill(mht_phi,wgt);
    h_mht_phi[this_category]->Fill(mht_phi,wgt);

    h_HT[0]->Fill(HT30,wgt);
    h_HT[this_category]->Fill(HT30,wgt);


    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      TLorentzVector myJet;
      myJet.SetPtEtaPhiE( eve->jet_pt_[iJet], eve->jet_eta_[iJet], eve->jet_phi_[iJet], eve->jet_energy_[iJet] );

      if( !firstJet ) sumJet = myJet;
      else            sumJet += myJet;

      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      double phi = eve->jet_phi_[iJet];

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      double dR = myJet.DeltaR(myLep);
      if( dR < 0.4 ) continue;

      double csv =  eve->jet_csv_[iJet];
      if( csv<0 && csv>-9 ) csv = -0.2;
      else if( csv < -5 )   csv = -0.4;

      if( csv > 1.0 ) csv = 1.0;

      h_jet_pt[0]->Fill(pt,wgt);
      h_jet_eta[0]->Fill(eta,wgt);
      h_jet_phi[0]->Fill(phi,wgt);
      h_jet_csv[0]->Fill(csv,wgt);

      h_jet_pt[this_category]->Fill(pt,wgt);
      h_jet_eta[this_category]->Fill(eta,wgt);
      h_jet_phi[this_category]->Fill(phi,wgt);
      h_jet_csv[this_category]->Fill(csv,wgt);
    }


  } // end loop over events

  std::cout << "**************************************************************" << std::endl;
  std::cout << "\t Number of raw events = " << numEvents_all << std::endl;
  std::cout << "\t Number of gen weighted events = " << numEvents_wgt_gen << std::endl;
  std::cout << "\t Number of lumi weighted events = " << numEvents_wgt_lumi << std::endl;
  std::cout << "\t Number of gen * lumi weighted events = " << numEvents_wgt << std::endl;
  std::cout << "**************************************************************" << std::endl;



  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();

}


double reweightPU( int nPU ){

  double PUscale[50];

  PUscale[0] = 6.83938;
  PUscale[1] = 2.30668;
  PUscale[2] = 1.91805;
  PUscale[3] = 1.83689;
  PUscale[4] = 1.89282;
  PUscale[5] = 1.94686;
  PUscale[6] = 1.97328;
  PUscale[7] = 1.93306;
  PUscale[8] = 1.82852;
  PUscale[9] = 1.67631;
  PUscale[10] = 1.48487;
  PUscale[11] = 1.27255;
  PUscale[12] = 1.0549;
  PUscale[13] = 0.845621;
  PUscale[14] = 0.665044;
  PUscale[15] = 0.508729;
  PUscale[16] = 0.389761;
  PUscale[17] = 0.289674;
  PUscale[18] = 0.217065;
  PUscale[19] = 0.161309;
  PUscale[20] = 0.118954;
  PUscale[21] = 0.0885512;
  PUscale[22] = 0.0675094;
  PUscale[23] = 0.0491225;
  PUscale[24] = 0.0386613;
  PUscale[25] = 0.0296305;
  PUscale[26] = 0.0252724;
  PUscale[27] = 0.0175392;
  PUscale[28] = 0.0135296;
  PUscale[29] = 0.011069;
  PUscale[30] = 0.00914227;
  PUscale[31] = 0.00730368;
  PUscale[32] = 0.00265843;
  PUscale[33] = 0.0021888;
  PUscale[34] = 0;
  PUscale[35] = 0;
  PUscale[36] = 0;
  PUscale[37] = 0;
  PUscale[38] = 0;
  PUscale[39] = 0;
  PUscale[40] = 0;
  PUscale[41] = 0;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;

  if( nPU>49 ) nPU = 49;

  return PUscale[nPU];
}


float DeltaR(float eta1,float phi1,float eta2,float phi2){
  float deltaPhi = TMath::Abs(phi1-phi2);
  float deltaEta = eta1-eta2;
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


/// algos / supporting functions:

void fillCSVhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}

double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  // force using no systematic here
  //iSys = 0;

  int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPt.size()); iJet++ ){

    double csv = jetCSV[iJet];
    double pt = jetPt[iJet];
    double jetAbsEta = fabs(jetEta[iJet]);
    int flavor = jetFlavor[iJet];

    int iPt = -1; int iEta = -1;
    if (pt >=19.99 && pt<30) iPt = 0;
    else if (pt >=30 && pt<40) iPt = 1;
    else if (pt >=40 && pt<60) iPt = 2;
    else if (pt >=60 && pt<100) iPt = 3;
    else if (pt >=100)          iPt = 4;
    //else if (pt >=100 && pt<160) iPt = 4;
    //else if (pt >=160 && pt<10000) iPt = 5;

    if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;

    //std::cout << " flavor = " << flavor << ", csv = " << csv << ", iPt = " << iPt << ", iEta = " << iEta << ", iSysHF = " << iSysHF << ", iSysC = " << iSysC << ", iSysLF = " << iSysLF << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, csv, iCSVWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			    flavor, pt, jetAbsEta, csv, iCSVWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
      // if( iSysLF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}
