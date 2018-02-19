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

#include "TriggerStudyEventVars.h"
//#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"
//#include "TriggerRun2/TriggerAnalyzer/interface/BTagCalibrationStandalone.h"

// #include "CondFormats/BTauObjects/interface/BTagCalibration.h"
// #include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;

double reweightPU( int nPU, int iSys );
float DeltaR(float eta1,float phi1,float eta2,float phi2);

// ------------ csv applying functions -------------
void fillCSVhistos(TFile *fileHF, TFile *fileLF);
void fillCMVAhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
//double get_csv_wgt_old( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
double get_cmva_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCMVA, vint jetFlavor, int iSys, double &cmvaWgtHF, double &cmvaWgtLF, double &cmvaWgtCF );

// double get_btv_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, TString sys_name, bool use_csv, bool verbose );

// void fillCSVEffhistos(TFile *file);
// double get_csv_efficiency( double jetPt, double jetEta, int jetFlavor );


int PtBinsHF_ = 5;


// CSV reweighting
TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];

// old CSV reweighting, for comparison
// TH1D* h_csv_wgt_hf_old[9][5];
// TH1D* c_csv_wgt_hf_old[9][5];
// TH1D* h_csv_wgt_lf_old[9][4][3];

// CSV reweighting
TH1D* h_cmva_wgt_hf[9][5];
TH1D* c_cmva_wgt_hf[9][5];
TH1D* h_cmva_wgt_lf[9][4][3];

// BTV efficiency
// TH2D* h_a_jet_pt_eta_eff_ = NULL;
// TH2D* h_b_jet_pt_eta_eff_ = NULL;
// TH2D* h_c_jet_pt_eta_eff_ = NULL;
// TH2D* h_l_jet_pt_eta_eff_ = NULL;

// TH2D* h_a_jet_pt_eta_all_ = NULL;
// TH2D* h_b_jet_pt_eta_all_ = NULL;
// TH2D* h_c_jet_pt_eta_all_ = NULL;
// TH2D* h_l_jet_pt_eta_all_ = NULL;


// setup calibration readers



//*****************************************************************************

void csv_normFix_treeReader( int insample=2500, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, int ttCat_=-1, bool useHTbins_=false, bool useCondor_=false ) {

  ////////////
  // b-tagging WPs --- keep these up-to-date
  ///////////
  double csvWPL = 0.1522; //0.5426;
  double csvWPM = 0.4941; //0.8484;
  double csvWPT = 0.8001; //0.9535;
  
  double cMVAWPL = 0.5803;//-0.5884;
  double cMVAWPM = 0.8838;//0.4432;
  double cMVAWPT = 0.9693;//0.9432;

  //--------
  bool rmPUJet = true;

  ///////////
  //////////
   TString inputFileHF = "data/csv_rwt_fit_hf_v2_final_2018_2_12test.root";
   TString inputFileLF = "data/csv_rwt_fit_lf_v2_final_2018_2_12test.root";

   if( useCondor_ ){
    inputFileHF = "csv_rwt_fit_hf_v2_final_2018_2_12test.root";
    inputFileLF = "csv_rwt_fit_lf_v2_final_2018_2_12test.root";
   }
   TFile* f_CSVwgt_HF = new TFile(inputFileHF);
   TFile* f_CSVwgt_LF = new TFile(inputFileLF);

  
  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

  // TFile* f_CSVwgt_HF_old = new TFile("data/csv_rwt_fit_hf_2015_11_20.root");
  // TFile* f_CSVwgt_LF_old = new TFile("data/csv_rwt_fit_lf_2015_11_20.root");

  // fillCSVhistos_old(f_CSVwgt_HF_old, f_CSVwgt_LF_old);

   TString inputFileHFCMVA = "data/cmva_rwt_fit_hf_v0_final_2018_2_13.root";
   TString inputFileLFCMVA = "data/cmva_rwt_fit_lf_v0_final_2018_2_13.root";

   if( useCondor_ ){
    inputFileHFCMVA = "cmva_rwt_fit_hf_v0_final_2018_2_13.root";
    inputFileLFCMVA = "cmva_rwt_fit_lf_v0_final_2018_2_13.root";
   }

  TFile* f_CMVAwgt_HF = new TFile(inputFileHFCMVA);
  TFile* f_CMVAwgt_LF = new TFile(inputFileLFCMVA);


  
  fillCMVAhistos(f_CMVAwgt_HF, f_CMVAwgt_LF);


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
    mySample_nGen_ = 97994442;//116591749;//25357774;//25446993;
    mySample_sampleName_ = "ttbar_powheg";
    // ext
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/lwming/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/10thFeb_trigger_csvRWT_13TeV/180210_151751/0000/");
    //    mySample_inputDirs_.push_back("/afs/cern.ch/work/l/lwming/csvRWT13TeV/");

    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_powheg";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_powheg";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_powheg";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_powheg";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_powheg";
    }
  }


  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/ttHbb_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());
  if( useCondor_ ) histofilename = Form("ttHbb_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());


  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    //if( insample!=2500 ) chain->Add(treefilename.c_str());
    chain->Add(treefilename.c_str());
  }


  
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

  TH1::SetDefaultSumw2();

  bool verbose_ = !true;

  //////////////////////////////////////////////////////////////////////////
  ///  Histograms
  //////////////////////////////////////////////////////////////////////////

  // int MaxNjet = 12;
  // int MaxNbtag = 5;

  
  // TH1D* h_numTruePVs = new TH1D("h_numTruePVs",";Number of True PVs", 50, 0, 50 );

  // TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );


  // std::vector<TString> cat_labels;
  // cat_labels.push_back("incl4j2t");
  // cat_labels.push_back("4j2t");
  // cat_labels.push_back("5j2t");
  // cat_labels.push_back("6j2t");
  // cat_labels.push_back("4j3t");
  // cat_labels.push_back("5j3t");
  // cat_labels.push_back("6j3t");
  // cat_labels.push_back("4j4t");
  // cat_labels.push_back("5j4t");
  // cat_labels.push_back("6j4t");

  // int NumCat = int(cat_labels.size());

  std::vector<TString> sys_cat_labels;
  sys_cat_labels.push_back("");               //0
  sys_cat_labels.push_back("_lepIdSFUp");     //1
  sys_cat_labels.push_back("_lepIdSFDown");   //2
  sys_cat_labels.push_back("_PUUp");          //3
  sys_cat_labels.push_back("_PUDown");        //4
  sys_cat_labels.push_back("_JERUp");         //5
  sys_cat_labels.push_back("_JERDown");       //6
  sys_cat_labels.push_back("_JESUp");         //7
  sys_cat_labels.push_back("_JESDown");       //8
  sys_cat_labels.push_back("_CSVLFUp");       //9
  sys_cat_labels.push_back("_CSVLFDown");     //10
  sys_cat_labels.push_back("_CSVHFUp");       //11
  sys_cat_labels.push_back("_CSVHFDown");     //12
  sys_cat_labels.push_back("_CSVHFStats1Up");     //13
  sys_cat_labels.push_back("_CSVHFStats1Down");   //14
  sys_cat_labels.push_back("_CSVHFStats2Up");     //15
  sys_cat_labels.push_back("_CSVHFStats2Down");   //16
  sys_cat_labels.push_back("_CSVLFStats1Up");     //17
  sys_cat_labels.push_back("_CSVLFStats1Down");   //18
  sys_cat_labels.push_back("_CSVLFStats2Up");     //19
  sys_cat_labels.push_back("_CSVLFStats2Down");   //20
  sys_cat_labels.push_back("_CSVCFErr1Up");     //21
  sys_cat_labels.push_back("_CSVCFErr1Down");   //22
  sys_cat_labels.push_back("_CSVCFErr2Up");     //23
  sys_cat_labels.push_back("_CSVCFErr2Down");   //24
  // sys_cat_labels.push_back("_muFUp");           //25
  // sys_cat_labels.push_back("_muFDown");         //26
  // sys_cat_labels.push_back("_muRUp");           //27
  // sys_cat_labels.push_back("_muRDown");         //28
  // sys_cat_labels.push_back("_muRmuFUp");        //29
  // sys_cat_labels.push_back("_muRmuFDown");      //30
  // sys_cat_labels.push_back("_NNPDFUp");        //31
  // sys_cat_labels.push_back("_NNPDFDown");      //32


  int NumSysCat = int(sys_cat_labels.size());


  TH1D* h_numEvents_perSys = new TH1D("h_numEvents_perSys",";Systematic", NumSysCat, 0-0.5, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV = new TH1D("h_numEvents_perSys_wgtCSV",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  // TH1D* h_numEvents_perSys_wgtCSV_4j2t = new TH1D("h_numEvents_perSys_wgtCSV_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  // TH1D* h_numEvents_perSys_wgtCSV_6j4t = new TH1D("h_numEvents_perSys_wgtCSV_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );


  TH1D* h_numEvents_perSys_wgtCMVA = new TH1D("h_numEvents_perSys_wgtCMVA",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  // TH1D* h_numEvents_perSys_wgtCMVA_4j2t = new TH1D("h_numEvents_perSys_wgtCMVA_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  // TH1D* h_numEvents_perSys_wgtCMVA_6j4t = new TH1D("h_numEvents_perSys_wgtCMVA_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );


  TH2D* h_hf_wgtCSV_perSys = new TH2D("h_hf_wgtCSV_perSys",";Systematic;hf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_lf_wgtCSV_perSys = new TH2D("h_lf_wgtCSV_perSys",";Systematic;lf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_cf_wgtCSV_perSys = new TH2D("h_cf_wgtCSV_perSys",";Systematic;cf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );

  
  int NumCuts = 5;

  /// Dilepton - Z control region


  ///// leq2j

  //// met30

  //// new


  //// met30


  // EleMu CR


  //// new



  //// Histograming 
  double metmax   = 500.;
  double lepPtMax = 300.;
  double jetptmax = 500.;
  double htmax    = 2000.;
  int NcsvBins = 107;//132;

  int NmetBins   = int( metmax/20. + 0.0001 );
  int NlepPtBins = int( lepPtMax/10. + 0.0001 );
  int NjetptBins = int( jetptmax/10. + 0.0001 );
  int NhtBins    = int( htmax/50. + 0.0001 );


  double HTmax = 1000.;
  int numHTbins = 1000;
  double L1HTmax = HTmax;
  int numL1HTTbins = numHTbins;



  std::vector<int> csv_systematics;
  csv_systematics.push_back(0);
  csv_systematics.push_back(7);
  csv_systematics.push_back(8);
  csv_systematics.push_back(9);
  csv_systematics.push_back(10);
  csv_systematics.push_back(11);
  csv_systematics.push_back(12);
  csv_systematics.push_back(13);
  csv_systematics.push_back(14);
  csv_systematics.push_back(15);
  csv_systematics.push_back(16);
  csv_systematics.push_back(17);
  csv_systematics.push_back(18);
  csv_systematics.push_back(19);
  csv_systematics.push_back(20);
  csv_systematics.push_back(21);
  csv_systematics.push_back(22);
  csv_systematics.push_back(23);
  csv_systematics.push_back(24);

  int numCSVSys = int(csv_systematics.size());


  TH1D* h_numEvents_Sys = new TH1D("h_numEvents_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );
  // TH1D* h_numEvents_4j2t_Sys = new TH1D("h_numEvents_4j2t_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );
  // TH1D* h_numEvents_4j_Sys = new TH1D("h_numEvents_4j_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );

  // TH1D* h_numEvents_Sys_PU = new TH1D("h_numEvents_Sys_PU",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );

  
  TH1D* h_numEvents_jet30_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet30_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCMVA_WP[numCSVSys];

  //
  TH1D* h_numEvents_jet30to50_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet30to50_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet50to70_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet50to70_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet70to100_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet70to100_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet100toInf_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet100toInf_nowgtCSV_WP[numCSVSys];

  //
  TH1D* h_numEvents_jet30to50_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30to50_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet50to70_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet50to70_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet70to100_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet70to100_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet100toInf_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet100toInf_nowgtCMVA_WP[numCSVSys];

  //
  TH1D* h_numEvents_lf_jet30_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_lf_jet30_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_lf_jet20_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_lf_jet20_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_lf_jet30_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_lf_jet30_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_lf_jet20_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_lf_jet20_nowgtCMVA_WP[numCSVSys];

  TH1D* h_wgt_csv_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_Sys[numCSVSys];

  // TH1D* h_wgt_csv_4j2t_Sys[numCSVSys];
  // TH1D* h_wgt_csv_hf_4j2t_Sys[numCSVSys];
  // TH1D* h_wgt_csv_lf_4j2t_Sys[numCSVSys];
  // TH1D* h_wgt_csv_cf_4j2t_Sys[numCSVSys];

  // TH1D* h_wgt_csv_4j_Sys[numCSVSys];
  // TH1D* h_wgt_csv_hf_4j_Sys[numCSVSys];
  // TH1D* h_wgt_csv_lf_4j_Sys[numCSVSys];
  // TH1D* h_wgt_csv_cf_4j_Sys[numCSVSys];


  ////
  TH1D* h_perjet_wgt_csv_Sys[5][3][numCSVSys];
  TH1D* h_perjet_wgt_csv_hf_Sys[5][1][numCSVSys];
  TH1D* h_perjet_wgt_csv_lf_Sys[4][3][numCSVSys];
  TH1D* h_perjet_wgt_csv_cf_Sys[5][1][numCSVSys];

  TH1D* h_perjet_wgt_cmva_Sys[5][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_hf_Sys[5][1][numCSVSys];
  TH1D* h_perjet_wgt_cmva_lf_Sys[4][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_cf_Sys[5][1][numCSVSys];

  for( int iSys=0; iSys<numCSVSys; iSys++ ){
    int useSys = csv_systematics[iSys];

    h_numEvents_jet30_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet20_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet20_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet20_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet20_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet30_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet20_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet20_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet20_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet20_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    //
    h_numEvents_jet30to50_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30to50_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30to50_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30to50_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet50to70_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet50to70_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet50to70_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet50to70_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet70to100_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet70to100_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet70to100_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet70to100_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet100toInf_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet100toInf_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet100toInf_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet100toInf_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    //
    h_numEvents_jet30to50_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30to50_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30to50_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30to50_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet50to70_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet50to70_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet50to70_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet50to70_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet70to100_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet70to100_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet70to100_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet70to100_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet100toInf_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet100toInf_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet100toInf_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet100toInf_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_wgt_csv_Sys[iSys] = new TH1D(Form("h_wgt_csv_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_Sys_%d",useSys),";weight", 300, 0., 3. );

    //
    h_numEvents_lf_jet30_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet30_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_lf_jet30_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet30_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_lf_jet20_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet20_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_lf_jet20_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet20_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_lf_jet30_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet30_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_lf_jet30_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet30_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_lf_jet20_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet20_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_lf_jet20_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_lf_jet20_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    // h_wgt_csv_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_hf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_lf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_cf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );

    // h_wgt_csv_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_hf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_lf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    // h_wgt_csv_cf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );


    
    for( int iPt=0; iPt < 5; iPt++ ){
      for( int iEta=0; iEta < 3; iEta++ ){
	h_perjet_wgt_csv_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_csv_hf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_hf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_csv_cf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_cf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iPt<4  ) h_perjet_wgt_csv_lf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_lf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );

	h_perjet_wgt_cmva_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_cmva_hf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_hf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_cmva_cf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_cf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iPt<4  ) h_perjet_wgt_cmva_lf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_lf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
      }
    }
  }



  //
  // Lepton plots

  // 
  // Jet plots
  //


  //
  // Energy sum plots
  //

  // std::vector<double> jet_pt_xbins;
  // for( int iBin=0; iBin<200; iBin++ ) jet_pt_xbins.push_back(double(iBin));
  // for( int iBin=200; iBin<400; iBin+=2 ) jet_pt_xbins.push_back(double(iBin));
  // for( int iBin=400; iBin<500; iBin+=10 ) jet_pt_xbins.push_back(double(iBin));
  // jet_pt_xbins.push_back(500);
  // int Njetbins = int(jet_pt_xbins.size());

  // double xbins_jetpt[311];
  // for( int iBin=0; iBin<Njetbins; iBin++ ) xbins_jetpt[iBin] = jet_pt_xbins[iBin];

  // TH2D* h_a_jet_pt_eta_csvM = new TH2D("h_a_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_b_jet_pt_eta_csvM = new TH2D("h_b_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_c_jet_pt_eta_csvM = new TH2D("h_c_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_l_jet_pt_eta_csvM = new TH2D("h_l_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );

  // TH2D* h_a_jet_pt_eta_all = new TH2D("h_a_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_b_jet_pt_eta_all = new TH2D("h_b_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_c_jet_pt_eta_all = new TH2D("h_c_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  // TH2D* h_l_jet_pt_eta_all = new TH2D("h_l_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );

  TH1D* h_jet_csv_b = new TH1D("h_jet_csv_b",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_c = new TH1D("h_jet_csv_c",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_l = new TH1D("h_jet_csv_l",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_a = new TH1D("h_jet_csv_a",";jet CSV", NcsvBins, -0.06, 1.01 );

  TH1D* h_jet_cmva_b = new TH1D("h_jet_cmva_b",";jet CMVA", 202, -1.01, 1.01 );
  TH1D* h_jet_cmva_c = new TH1D("h_jet_cmva_c",";jet CMVA", 202, -1.01, 1.01 );
  TH1D* h_jet_cmva_l = new TH1D("h_jet_cmva_l",";jet CMVA", 202, -1.01, 1.01 );
  TH1D* h_jet_cmva_a = new TH1D("h_jet_cmva_a",";jet CMVA", 202, -1.01, 1.01 );

  TH2D* h_b_csv_wgt = new TH2D("h_b_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );

  TH2D* h_c_csv_wgt = new TH2D("h_c_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );

  TH2D* h_l_csv_wgt = new TH2D("h_l_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );



  int NumHists = 6;
  int NumSysC = 5;
  TH1D* h_c_jet_csv[NumHists];
  TH1D* h_c_jet_csv_CharmCSVSF[NumSysC][NumHists];

  TH1D* h_c_jet_cmva[NumHists];
  TH1D* h_c_jet_cmva_CharmCMVASF[NumSysC][NumHists];

  int NfullcsvBins = 1000;
  for( int iPt=0; iPt<NumHists; iPt++ ){
    h_c_jet_csv[iPt] = new TH1D(Form("h_c_jet_csv_Pt%d_Eta%d", iPt, 0),";jet CSV", NfullcsvBins, -0.1, 1.1 );
    h_c_jet_cmva[iPt] = new TH1D(Form("h_c_jet_cmva_Pt%d_Eta%d", iPt, 0),";jet CMVA", NfullcsvBins, -1.1, 1.1 );
    for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
      h_c_jet_csv_CharmCSVSF[iSysC][iPt] = new TH1D(Form("h_c_jet_csv_CharmCSVSF_Pt%d_Eta%d_SysC%d", iPt, 0, iSysC),";jet CSV", NfullcsvBins, -0.1, 1.1 );
      h_c_jet_cmva_CharmCMVASF[iSysC][iPt] = new TH1D(Form("h_c_jet_cmva_CharmCMVASF_Pt%d_Eta%d_SysC%d", iPt, 0, iSysC),";jet CMVA", NfullcsvBins, -1.1, 1.1 );
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
  double numEvents_wgt_gen_lumi = 0;
  double numEvents_wgt_gen_lumi_pu = 0;
  double numEvents_wgt_gen_lumi_pu_csv = 0;

  double numEvents_wgt_gen_lumi_pu_4j2t = 0;
  double numEvents_wgt_gen_lumi_pu_csv_4j2t = 0;


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

    // std::cout << "     Event " << ievt << std::endl;

    chain->GetEntry(ievt);

    if( verbose_ ) std::cout << " ===> test 0 " << std::endl;
    // If using HT bins, only take low HT bins of inclusive samples
    // double lheHT = eve->lheHT_;
    // lheHT = std::min( lheHT, MaxLHEHT-0.00001 );
    
    // h_lheHT->Fill(lheHT);

    if( verbose_ ) std::cout << " ===> test 0.1 " << std::endl;

    // int additionalJetEventId = -99;
    // if( insample>=2500 && insample<=2504 ){
    //   additionalJetEventId = eve->additionalJetEventId_;

    //   bool keepTTbarEvent = true;
    //   if( ttCat_>=0 ){
    // 	if( ttCat_==0 && additionalJetEventId==0 ) keepTTbarEvent = true;
    // 	else if( ttCat_==1 && (additionalJetEventId>=41 && additionalJetEventId<=45) ) keepTTbarEvent = true;
    // 	else if( ttCat_==2 && additionalJetEventId==51 ) keepTTbarEvent = true;
    // 	else if( ttCat_==3 && additionalJetEventId==52 ) keepTTbarEvent = true;
    // 	else if( ttCat_==4 && (additionalJetEventId>=53 && additionalJetEventId<=55) ) keepTTbarEvent = true;
    // 	else keepTTbarEvent = false;
    //   }

    //   // if( keepTTbarEvent==false && ttCat_>=0 ){
    //   // 	std::cout << " ERROR!! keepTTbarEvent == false! additionalJetEventId = " << additionalJetEventId << " and ttCat_ = " << ttCat_ << std::endl;
    //   // }

    //   if( !keepTTbarEvent ) continue;
    //   //std::cout << " I have made it past the gate! " << std::endl;
    // }

    // h_additionalJetEventId->Fill(additionalJetEventId);

    // int numTruePVs = eve->numTruePVs_;

    // //double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numPVs);
    // double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,0);
    // double wgt_pu_up = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,1);
    // double wgt_pu_down = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,-1);
    if( verbose_ ) std::cout << " ===> test 0.2 " << std::endl;

    double wgt_pu = 1.;
    double wgt_pu_up = 1.;
    double wgt_pu_down = 1.;
    
    double wgt_gen = ( insample > 0 ) ? eve->wgt_generator_ : 1;
    wgt_gen = ( wgt_gen > 0 ) ? 1. : -1.;
    //if( ievt<3 ){
    // std::vector<double> lhe_event_weights = eve->LHEEvent_weights_;
    // //std::cout << "\t originalXWGTUP = " << eve->originalXWGTUP_ << std::endl;
    // for( unsigned int iwgt = 0; iwgt < lhe_event_weights.size(); iwgt++ ){
    //   h_numEvents_lheWeight->Fill(iwgt, lhe_event_weights[iwgt]/eve->originalXWGTUP_);
    // 	//printf(" weight %5d = %f,\t ratio = %4.5f\n", iwgt, lhe_event_weights[iwgt], lhe_event_weights[iwgt]/eve->originalXWGTUP_);
    // }
      //printf(" Event %3lld: additionalJetEventId = %3d \n", ievt, additionalJetEventId);
    //}



    //
    // Lepton selection
    //

    std::vector<int> ind_ele;
    std::vector<int> ind_mu;
    std::vector<int> ind_ele_loose;
    std::vector<int> ind_mu_loose;
    for( int iLep=0; iLep<int(eve->lepton_pt_.size()); iLep++ ){
      bool isMuon = ( eve->lepton_isMuon_[iLep]==1 );

      bool isTrigMVAM = ( eve->lepton_isTrigMVAM_[iLep]==1 );
      bool isTight = ( eve->lepton_isTight_[iLep]==1 );
      bool isLoose = ( eve->lepton_isLoose_[iLep]==1 );
      bool isCrack = ( eve->lepton_inCrack_[iLep]==1 );

      double pt  = eve->lepton_pt_[iLep];
      double eta = eve->lepton_eta_[iLep];

      double relIso = eve->lepton_relIso_[iLep];
      double relIsoR04 = eve->lepton_relIsoR04_[iLep];

      ////
      if( isMuon ){
	if( pt > 25 && abs(eta)<2.1 && isTight && relIsoR04 < 0.15 ) ind_mu.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isLoose && relIsoR04 < 0.25 ) ind_mu_loose.push_back(iLep);
      }

      if( !isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isTrigMVAM && !isCrack && relIso < 0.15 ) ind_ele.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isTrigMVAM && !isCrack && relIso < 0.15 ) ind_ele_loose.push_back(iLep);
      }
    }

        if( verbose_ ) std::cout << " ===> test 0.3 " << std::endl;

    std::vector<double> rejetPts;
    std::vector<double> rejetEtas;
    std::vector<double> rejetCSVs;
    std::vector<int>    rejetFlavors;

    int njet = 0;
    int nbtag = 0;

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
       if( verbose_ ) std::cout << " ===> test 0.3.0.1 " << std::endl;
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      int flavor = eve->jet_hadronFlavour_[iJet];
      double csv =  eve->jet_csv_[iJet];
      double cmva = eve->jet_cmva_[iJet];
       if( verbose_ ) std::cout << " ===> test 0.3.0.2 " << std::endl;

      double jetAbsEta = fabs(eta);
      
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;
      if( rmPUJet && (eve->jet_PUID_passWPLoose_[iJet] != 1) ) continue;
      if( !(pt>20. && fabs(eta)<2.5) ) continue;
             if( verbose_ ) std::cout << " ===> test 0.3.0.3 " << std::endl;

      int iPt = -1, iEta = -1;
      if (pt >=19.99 && pt<30) iPt = 0;
      else if (pt >=30 && pt<40) iPt = 1;
      else if (pt >=40 && pt<60) iPt = 2;
      else if (pt >=60 && pt<100) iPt = 3;
      else if (pt >=100)          iPt = 3;
      //else if (pt >=100 && pt<160) iPt = 4;
      //else if (pt >=160 && pt<10000) iPt = 5;

      ////
      if(abs(flavor) == 5 || abs(flavor) == 4){
	if (pt >=19.99 && pt<30) iPt = 0;
	else if (pt >=30 && pt<50) iPt = 1;
	else if (pt >=50 && pt<70) iPt = 2;
	else if (pt >=70 && pt<100) iPt = 3;
	else if (pt >=100 && pt<160) iPt = 4;
	else if (pt >=160 && pt<10000) iPt = 4; //4
      }

      if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;
       if( verbose_ ) std::cout << " ===> test 0.3.0.4 " << std::endl;

      
      if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
      else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
      else if ( jetAbsEta>=1.6 && jetAbsEta<2.51 ) iEta = 2;
        if( verbose_ ) std::cout << " ===> test 0.3.1 " << std::endl;


      if( abs(flavor)==4 && fabs(eta)<2.5 && pt>20. ){
	h_c_jet_csv[iPt]->Fill(csv);
	h_c_jet_cmva[iPt]->Fill(cmva);
	for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
	  int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	  double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	  int useCMVABin = (cmva>=0.) ? c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva) : 1;
	  double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);

	  h_c_jet_csv_CharmCSVSF[iSysC][iPt]->Fill(csv,iCSVWgtC);
	  h_c_jet_cmva_CharmCMVASF[iSysC][iPt]->Fill(cmva,iCMVAWgtC);
	}
      }
       if( verbose_ ) std::cout << " ===> test 0.3.2 " << std::endl;

      if( !(pt>30. && fabs(eta)<2.5) ) continue;

      h_jet_csv_a->Fill(csv);
      if( abs(flavor)==5 )      h_jet_csv_b->Fill(csv);
      else if( abs(flavor)==4 ) h_jet_csv_c->Fill(csv);
      else                      h_jet_csv_l->Fill(csv);

      h_jet_cmva_a->Fill(cmva);
      if( abs(flavor)==5 )      h_jet_cmva_b->Fill(cmva);
      else if( abs(flavor)==4 ) h_jet_cmva_c->Fill(cmva);
      else                      h_jet_cmva_l->Fill(cmva);

      njet++;
      if( csv>csvWPM ) nbtag++;

      rejetPts.push_back(pt);
      rejetEtas.push_back(eta);
      rejetCSVs.push_back(csv);
      rejetFlavors.push_back(flavor);


      pt = std::min(pt, jetptmax-0.001);
      eta = fabs(eta);
      
      // h_a_jet_pt_eta_all->Fill(pt, eta);
      // if( abs(flavor)==5 )      h_b_jet_pt_eta_all->Fill(pt, eta);
      // else if( abs(flavor)==4 ) h_c_jet_pt_eta_all->Fill(pt, eta);
      // else                      h_l_jet_pt_eta_all->Fill(pt, eta);

      // if( csv>0.80 ){
      // 	h_a_jet_pt_eta_csvM->Fill(pt, eta);
      // 	if( abs(flavor)==5 )      h_b_jet_pt_eta_csvM->Fill(pt, eta);
      // 	else if( abs(flavor)==4 ) h_c_jet_pt_eta_csvM->Fill(pt, eta);
      // 	else                      h_l_jet_pt_eta_csvM->Fill(pt, eta);
      // }
      
    }
    if( verbose_ ) std::cout << " ===> test 0.4 " << std::endl;






    // JESup
    std::vector<double> rejetPts_JESup;
    std::vector<double> rejetEtas_JESup;
    std::vector<double> rejetCSVs_JESup;
    std::vector<int>    rejetFlavors_JESup;

    int njet_JESup = 0;
    int nbtag_JESup = 0;

    for( int iJet=0; iJet<int(eve->jet_JESup_pt_.size()); iJet++ ){
      double pt  = eve->jet_JESup_pt_[iJet];
      double eta = eve->jet_JESup_eta_[iJet];
      int flavor = eve->jet_JESup_hadronFlavour_[iJet];
      double csv =  eve->jet_JESup_csv_[iJet];

      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;
      if( rmPUJet && (eve->jet_JESup_PUID_passWPLoose_[iJet] != 1) ) continue;
      if( !(pt>30. && fabs(eta)<2.5) ) continue;

      njet_JESup++;
      if( csv>csvWPM ) nbtag_JESup++;

      rejetPts_JESup.push_back(pt);
      rejetEtas_JESup.push_back(eta);
      rejetCSVs_JESup.push_back(csv);
      rejetFlavors_JESup.push_back(flavor);
    }



    // JESdown
    std::vector<double> rejetPts_JESdown;
    std::vector<double> rejetEtas_JESdown;
    std::vector<double> rejetCSVs_JESdown;
    std::vector<int>    rejetFlavors_JESdown;

    int njet_JESdown = 0;
    int nbtag_JESdown = 0;

    for( int iJet=0; iJet<int(eve->jet_JESdown_pt_.size()); iJet++ ){
      double pt  = eve->jet_JESdown_pt_[iJet];
      double eta = eve->jet_JESdown_eta_[iJet];
      int flavor = eve->jet_JESdown_hadronFlavour_[iJet];
      double csv =  eve->jet_JESdown_csv_[iJet];
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;
      if( rmPUJet && (eve->jet_JESdown_PUID_passWPLoose_[iJet] != 1) ) continue;
      if( !(pt>30. && fabs(eta)<2.5) ) continue;

      njet_JESdown++;
      if( csv>csvWPM ) nbtag_JESdown++;

      rejetPts_JESdown.push_back(pt);
      rejetEtas_JESdown.push_back(eta);
      rejetCSVs_JESdown.push_back(csv);
      rejetFlavors_JESdown.push_back(flavor);
    }

    if( verbose_ ) std::cout << " ===> test 1 " << std::endl;

    int numPVs = eve->numPVs_;

    double wgt_csv_noSys = 1;
    double wgt_csv_hf_noSys = 1;
    double wgt_csv_lf_noSys = 1;
    double wgt_csv_cf_noSys = 1;
    if( verbose_ ) std::cout << " ===> test 0.6 " << std::endl;

    for( int useSys=0; useSys<int(csv_systematics.size()); useSys++ ){
      int mySys = csv_systematics[useSys];

      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;

      double temp_wgt_csv = 1.;//( insample<0 ) ? 1 : get_csv_wgt(rejetPts, rejetEtas, rejetCSVs, rejetFlavors, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      if( !(insample < 0) ){
      	if( mySys==7 )      temp_wgt_csv = get_csv_wgt(rejetPts_JESup, rejetEtas_JESup, rejetCSVs_JESup, rejetFlavors_JESup, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      	else if( mySys==8 ) temp_wgt_csv = get_csv_wgt(rejetPts_JESdown, rejetEtas_JESdown, rejetCSVs_JESdown, rejetFlavors_JESdown, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      	else                temp_wgt_csv = get_csv_wgt(rejetPts, rejetEtas, rejetCSVs, rejetFlavors, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      }

      if( useSys==0 ){
	wgt_csv_noSys = temp_wgt_csv;
	wgt_csv_hf_noSys = wgt_csv_hf;
	wgt_csv_lf_noSys = wgt_csv_lf;
	wgt_csv_cf_noSys = wgt_csv_cf;


	// h_numJet_nolepreq->Fill(njet,wgt_gen*wgt_pu);
	// h_numJet_nolepreq_wgtCSV->Fill(njet,wgt_gen*wgt_pu*temp_wgt_csv);
      }

      
      h_numEvents_Sys->Fill(0.5+mySys,temp_wgt_csv);
      // h_numEvents_Sys_PU->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);

      h_wgt_csv_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
      h_wgt_csv_hf_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
      h_wgt_csv_lf_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
      h_wgt_csv_cf_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);


      int temp_njet = njet;
      int temp_nbtag = nbtag;

      if( insample < 0 ){
	if( mySys==7 ){
	  temp_njet = njet_JESup;
	  temp_nbtag = nbtag_JESup;
	}
	else if( mySys==8 ){
	  temp_njet = njet_JESdown;
	  temp_nbtag = nbtag_JESdown;
	}
	else {
	  temp_njet = njet;
	  temp_nbtag = nbtag;
	}
      }

      // if( temp_njet>=4 && temp_nbtag>=2 ){
      // 	h_numEvents_4j2t_Sys->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);
      // 	h_wgt_csv_4j2t_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
      // 	h_wgt_csv_hf_4j2t_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
      // 	h_wgt_csv_lf_4j2t_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
      // 	h_wgt_csv_cf_4j2t_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);
      // }
      
      // if( temp_njet>=4 ){
      // 	h_numEvents_4j_Sys->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);
      // 	h_wgt_csv_4j_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
      // 	h_wgt_csv_hf_4j_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
      // 	h_wgt_csv_lf_4j_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
      // 	h_wgt_csv_cf_4j_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);

      // }



      //// NEW

      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_cmva = eve->jet_cmva_;
      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;
      vint use_jet_PUID_passWPLoose = eve->jet_PUID_passWPLoose_;
      if( mySys==7 ){
	use_jet_pt = eve->jet_JESup_pt_;
	use_jet_eta = eve->jet_JESup_eta_;
	use_jet_csv = eve->jet_JESup_csv_;
	use_jet_cmva = eve->jet_JESup_cmva_;
	use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;
	use_jet_PUID_passWPLoose = eve->jet_JESup_PUID_passWPLoose_;
      }
      else if( mySys==8 ){
	use_jet_pt = eve->jet_JESdown_pt_;
	use_jet_eta = eve->jet_JESdown_eta_;
	use_jet_csv = eve->jet_JESdown_csv_;
	use_jet_cmva = eve->jet_JESdown_cmva_;
	use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;
	use_jet_PUID_passWPLoose = eve->jet_JESdown_PUID_passWPLoose_;
      }
      else {
	use_jet_pt = eve->jet_pt_;
	use_jet_eta = eve->jet_eta_;
	use_jet_csv = eve->jet_csv_;
	use_jet_cmva = eve->jet_cmva_;
	use_jet_hadronFlavour = eve->jet_hadronFlavour_;
	use_jet_PUID_passWPLoose = eve->jet_PUID_passWPLoose_;
      }


      double wgt_jet30_CSVL = 0.;
      double wgt_jet30_CSVM = 0.;
      double wgt_jet30_CSVT = 0.;

      double nowgt_jet30_CSVL = 0.;
      double nowgt_jet30_CSVM = 0.;
      double nowgt_jet30_CSVT = 0.;

      double wgt_jet20_CSVL = 0.;
      double wgt_jet20_CSVM = 0.;
      double wgt_jet20_CSVT = 0.;

      double nowgt_jet20_CSVL = 0.;
      double nowgt_jet20_CSVM = 0.;
      double nowgt_jet20_CSVT = 0.;


      double wgt_jet30_CMVAL = 0.;
      double wgt_jet30_CMVAM = 0.;
      double wgt_jet30_CMVAT = 0.;

      double nowgt_jet30_CMVAL = 0.;
      double nowgt_jet30_CMVAM = 0.;
      double nowgt_jet30_CMVAT = 0.;

      double wgt_jet20_CMVAL = 0.;
      double wgt_jet20_CMVAM = 0.;
      double wgt_jet20_CMVAT = 0.;

      double nowgt_jet20_CMVAL = 0.;
      double nowgt_jet20_CMVAM = 0.;
      double nowgt_jet20_CMVAT = 0.;

      ///////
      double wgt_jet30to50_CSVL = 0.;
      double wgt_jet30to50_CSVM = 0.;
      double wgt_jet30to50_CSVT = 0.;

      double nowgt_jet30to50_CSVL = 0.;
      double nowgt_jet30to50_CSVM = 0.;
      double nowgt_jet30to50_CSVT = 0.;


      double wgt_jet50to70_CSVL = 0.;
      double wgt_jet50to70_CSVM = 0.;
      double wgt_jet50to70_CSVT = 0.;

      double nowgt_jet50to70_CSVL = 0.;
      double nowgt_jet50to70_CSVM = 0.;
      double nowgt_jet50to70_CSVT = 0.;

      double wgt_jet70to100_CSVL = 0.;
      double wgt_jet70to100_CSVM = 0.;
      double wgt_jet70to100_CSVT = 0.;

      double nowgt_jet70to100_CSVL = 0.;
      double nowgt_jet70to100_CSVM = 0.;
      double nowgt_jet70to100_CSVT = 0.;


      double wgt_jet100toInf_CSVL = 0.;
      double wgt_jet100toInf_CSVM = 0.;
      double wgt_jet100toInf_CSVT = 0.;

      double nowgt_jet100toInf_CSVL = 0.;
      double nowgt_jet100toInf_CSVM = 0.;
      double nowgt_jet100toInf_CSVT = 0.;

      /////
      double wgt_jet30to50_CMVAL = 0.;
      double wgt_jet30to50_CMVAM = 0.;
      double wgt_jet30to50_CMVAT = 0.;

      double nowgt_jet30to50_CMVAL = 0.;
      double nowgt_jet30to50_CMVAM = 0.;
      double nowgt_jet30to50_CMVAT = 0.;


      double wgt_jet50to70_CMVAL = 0.;
      double wgt_jet50to70_CMVAM = 0.;
      double wgt_jet50to70_CMVAT = 0.;

      double nowgt_jet50to70_CMVAL = 0.;
      double nowgt_jet50to70_CMVAM = 0.;
      double nowgt_jet50to70_CMVAT = 0.;

      double wgt_jet70to100_CMVAL = 0.;
      double wgt_jet70to100_CMVAM = 0.;
      double wgt_jet70to100_CMVAT = 0.;

      double nowgt_jet70to100_CMVAL = 0.;
      double nowgt_jet70to100_CMVAM = 0.;
      double nowgt_jet70to100_CMVAT = 0.;


      double wgt_jet100toInf_CMVAL = 0.;
      double wgt_jet100toInf_CMVAM = 0.;
      double wgt_jet100toInf_CMVAT = 0.;

      double nowgt_jet100toInf_CMVAL = 0.;
      double nowgt_jet100toInf_CMVAM = 0.;
      double nowgt_jet100toInf_CMVAT = 0.;

      //light flavor
      double wgt_lf_jet30_CSVL = 0.;
      double wgt_lf_jet30_CSVM = 0.;
      double wgt_lf_jet30_CSVT = 0.;

      double nowgt_lf_jet30_CSVL = 0.;
      double nowgt_lf_jet30_CSVM = 0.;
      double nowgt_lf_jet30_CSVT = 0.;

      double wgt_lf_jet20_CSVL = 0.;
      double wgt_lf_jet20_CSVM = 0.;
      double wgt_lf_jet20_CSVT = 0.;

      double nowgt_lf_jet20_CSVL = 0.;
      double nowgt_lf_jet20_CSVM = 0.;
      double nowgt_lf_jet20_CSVT = 0.;


      double wgt_lf_jet30_CMVAL = 0.;
      double wgt_lf_jet30_CMVAM = 0.;
      double wgt_lf_jet30_CMVAT = 0.;

      double nowgt_lf_jet30_CMVAL = 0.;
      double nowgt_lf_jet30_CMVAM = 0.;
      double nowgt_lf_jet30_CMVAT = 0.;

      double wgt_lf_jet20_CMVAL = 0.;
      double wgt_lf_jet20_CMVAM = 0.;
      double wgt_lf_jet20_CMVAT = 0.;

      double nowgt_lf_jet20_CMVAL = 0.;
      double nowgt_lf_jet20_CMVAM = 0.;
      double nowgt_lf_jet20_CMVAT = 0.;


      //////
      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];
	int flavor = use_jet_hadronFlavour[iJet];
	double csv =  use_jet_csv[iJet];
	double cmva =  use_jet_cmva[iJet];
	int PUID_passWPLoose = use_jet_PUID_passWPLoose[iJet];
	double jetAbsEta = fabs(eta);
      
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;
	if( cmva < 0.0 ) cmva = -0.05;
	if( cmva > 1.0 ) cmva = 1.0;

	if( rmPUJet && (PUID_passWPLoose != 1) ) continue;
	if( !(pt>20. && fabs(eta)<2.5) ) continue;
      
	int iPt = -1, iEta = -1;
	if (pt >=19.99 && pt<30) iPt = 0;
	else if (pt >=30 && pt<40) iPt = 1;
	else if (pt >=40 && pt<60) iPt = 2;
	else if (pt >=60 && pt<100) iPt = 3;
	else if (pt >=100)          iPt = 3;
	//else if (pt >=100 && pt<160) iPt = 4;
	//else if (pt >=160 && pt<10000) iPt = 5;

	////
	if(abs(flavor) == 5 || abs(flavor) == 4){
	  if (pt >=19.99 && pt<30) iPt = 0;
	  else if (pt >=30 && pt<50) iPt = 1;
	  else if (pt >=50 && pt<70) iPt = 2;
	  else if (pt >=70 && pt<100) iPt = 3;
	  else if (pt >=100 && pt<160) iPt = 4;
	  else if (pt >=160 && pt<10000) iPt = 4; //4
	}

	if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

      
	if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
	else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
	else if ( jetAbsEta>=1.6 && jetAbsEta<2.51 ) iEta = 2;

	int iSysHF = 0;
	switch(mySys){
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
	switch(mySys){
	case 21: iSysC=1; break;
	case 22: iSysC=2; break;
	case 23: iSysC=3; break;
	case 24: iSysC=4; break;
	default : iSysC = 0; break;
	}

	int iSysLF = 0;
	switch(mySys){
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
  
	if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;
	if( verbose_ ) std::cout << " ===> 2 test 1.3 jet " << iJet << std::endl;

	if (abs(flavor) == 5 ){
	  if( iEta>=1 ) iEta=0;
	  int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
	  double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_hf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtHF,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtHF,wgt_gen);

	  if( verbose_ ) std::cout << " ===> 2 test 1.4 jet " << iJet << std::endl;
	  int useCMVABin = (cmva>=0.) ? h_cmva_wgt_hf[iSysHF][iPt]->FindBin(cmva) : 1;
	  double iCMVAWgtHF = h_cmva_wgt_hf[iSysHF][iPt]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.5 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_hf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtHF,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtHF,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.6 jet " << iJet << std::endl;

	  if( csv > csvWPL && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVL += iCSVWgtHF;
	    nowgt_jet30_CSVL += 1.0;
	  }
	  if( csv > csvWPM && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVM += iCSVWgtHF;
	    nowgt_jet30_CSVM += 1.0;
	  }
	  if( csv > csvWPT && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVT += iCSVWgtHF;
	    nowgt_jet30_CSVT += 1.0;
	  }

	  
	  if( csv > csvWPL && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVL += iCSVWgtHF;
	    nowgt_jet20_CSVL += 1.0;
	  }
	  if( csv > csvWPM && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVM += iCSVWgtHF;
	    nowgt_jet20_CSVM += 1.0;
	  }
	  if( csv > csvWPT && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVT += iCSVWgtHF;
	    nowgt_jet20_CSVT += 1.0;
	  }


	  //---------------------
	  if( pt>=30 && pt<50 && iCSVWgtHF>0. ){
	    if( csv > csvWPL ){
	      wgt_jet30to50_CSVL += iCSVWgtHF;
	      nowgt_jet30to50_CSVL += 1.0;
	    }
	    if( csv > csvWPM ){
	      wgt_jet30to50_CSVM += iCSVWgtHF;
	      nowgt_jet30to50_CSVM += 1.0;
	    }
	    if( csv > csvWPT ){
	      wgt_jet30to50_CSVT += iCSVWgtHF;
	      nowgt_jet30to50_CSVT += 1.0;
	    }
	  }

	  

	  if( pt>=50 && pt<70 && iCSVWgtHF>0. ){
	    if( csv > csvWPL ){
	      wgt_jet50to70_CSVL += iCSVWgtHF;
	      nowgt_jet50to70_CSVL += 1.0;
	    }
	    if( csv > csvWPM ){
	      wgt_jet50to70_CSVM += iCSVWgtHF;
	      nowgt_jet50to70_CSVM += 1.0;
	    }
	    if( csv > csvWPT ){
	      wgt_jet50to70_CSVT += iCSVWgtHF;
	      nowgt_jet50to70_CSVT += 1.0;
	    }
	  }

	  if( pt>=70 && pt<100 && iCSVWgtHF>0. ){
	    if( csv > csvWPL ){
	      wgt_jet70to100_CSVL += iCSVWgtHF;
	      nowgt_jet70to100_CSVL += 1.0;
	    }
	    if( csv > csvWPM ){
	      wgt_jet70to100_CSVM += iCSVWgtHF;
	      nowgt_jet70to100_CSVM += 1.0;
	    }
	    if( csv > csvWPT ){
	      wgt_jet70to100_CSVT += iCSVWgtHF;
	      nowgt_jet70to100_CSVT += 1.0;
	    }
	  }
 

	  if( pt>=100 && iCSVWgtHF>0. ){
	    if( csv > csvWPL ){
	      wgt_jet100toInf_CSVL += iCSVWgtHF;
	      nowgt_jet100toInf_CSVL += 1.0;
	    }
	    if( csv > csvWPM ){
	      wgt_jet100toInf_CSVM += iCSVWgtHF;
	      nowgt_jet100toInf_CSVM += 1.0;
	    }
	    if( csv > csvWPT ){
	      wgt_jet100toInf_CSVT += iCSVWgtHF;
	      nowgt_jet100toInf_CSVT += 1.0;
	    }
	  }

	  /////////////////////////
	  /////////////////////////

	  if( cmva > cMVAWPL && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAL += iCMVAWgtHF;
	    nowgt_jet30_CMVAL += 1.0;
	  }
	  if( cmva > cMVAWPM && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAM += iCMVAWgtHF;
	    nowgt_jet30_CMVAM += 1.0;
	  }
	  if( cmva > cMVAWPT && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAT += iCMVAWgtHF;
	    nowgt_jet30_CMVAT += 1.0;
	  }

	  
	  if( cmva > cMVAWPL && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAL += iCMVAWgtHF;
	    nowgt_jet20_CMVAL += 1.0;
	  }
	  if( cmva > cMVAWPM && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAM += iCMVAWgtHF;
	    nowgt_jet20_CMVAM += 1.0;
	  }
	  if( cmva > cMVAWPT && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAT += iCMVAWgtHF;
	    nowgt_jet20_CMVAT += 1.0;
	  }


	  
	  //---------------------
	  if( pt>=30 && pt<50 && iCMVAWgtHF>0. ){
	    if( cmva > cMVAWPL ){
	      wgt_jet30to50_CMVAL += iCMVAWgtHF;
	      nowgt_jet30to50_CMVAL += 1.0;
	    }
	    if( cmva > cMVAWPM ){
	      wgt_jet30to50_CMVAM += iCMVAWgtHF;
	      nowgt_jet30to50_CMVAM += 1.0;
	    }
	    if( cmva > cMVAWPT ){
	      wgt_jet30to50_CMVAT += iCMVAWgtHF;
	      nowgt_jet30to50_CMVAT += 1.0;
	    }
	  }

	  

	  if( pt>=50 && pt<70 && iCMVAWgtHF>0. ){
	    if( cmva > cMVAWPL ){
	      wgt_jet50to70_CMVAL += iCMVAWgtHF;
	      nowgt_jet50to70_CMVAL += 1.0;
	    }
	    if( cmva > cMVAWPM ){
	      wgt_jet50to70_CMVAM += iCMVAWgtHF;
	      nowgt_jet50to70_CMVAM += 1.0;
	    }
	    if( cmva > cMVAWPT ){
	      wgt_jet50to70_CMVAT += iCMVAWgtHF;
	      nowgt_jet50to70_CMVAT += 1.0;
	    }
	  }

	  if( pt>=70 && pt<100 && iCMVAWgtHF>0. ){
	    if( cmva > cMVAWPL ){
	      wgt_jet70to100_CMVAL += iCMVAWgtHF;
	      nowgt_jet70to100_CMVAL += 1.0;
	    }
	    if( cmva > cMVAWPM ){
	      wgt_jet70to100_CMVAM += iCMVAWgtHF;
	      nowgt_jet70to100_CMVAM += 1.0;
	    }
	    if( cmva > cMVAWPT ){
	      wgt_jet70to100_CMVAT += iCMVAWgtHF;
	      nowgt_jet70to100_CMVAT += 1.0;
	    }
	  }
 

	  if( pt>=100 && iCMVAWgtHF>0. ){
	    if( cmva > cMVAWPL ){
	      wgt_jet100toInf_CMVAL += iCMVAWgtHF;
	      nowgt_jet100toInf_CMVAL += 1.0;
	    }
	    if( cmva > cMVAWPM ){
	      wgt_jet100toInf_CMVAM += iCMVAWgtHF;
	      nowgt_jet100toInf_CMVAM += 1.0;
	    }
	    if( cmva > cMVAWPT ){
	      wgt_jet100toInf_CMVAT += iCMVAWgtHF;
	      nowgt_jet100toInf_CMVAT += 1.0;
	    }
	  }



	  
	}
	else if( abs(flavor) == 4 ){
	  if( iEta>=1 ) iEta=0;
	  int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	  double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_cf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtC,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtC,wgt_gen);

	  
	  if( verbose_ ) std::cout << " ===> 2 test 1.7 jet " << iJet << std::endl;
	  int useCMVABin = (cmva>=0.) ? c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva) : 1;
	  double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.8 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_cf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtC,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtC,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.9 jet " << iJet << std::endl;
	}
	else {
	  if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
	  int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
	  double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_lf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtLF,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtLF,wgt_gen);

	  
	  if( verbose_ ) std::cout << " ===> 2 test 1.10 jet " << iJet << std::endl;
	  int useCMVABin = (cmva>=0.) ? h_cmva_wgt_lf[iSysLF][iPt][iEta]->FindBin(cmva) : 1;
	  double iCMVAWgtLF = h_cmva_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.11 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_lf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtLF,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtLF,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.12 jet " << iJet << std::endl;


	  if( csv < csvWPL && pt > 30 && iCSVWgtLF>0. ){
	    wgt_lf_jet30_CSVL += iCSVWgtLF;
	    nowgt_lf_jet30_CSVL += 1.0;
	  }
	  if( csv < csvWPM && pt > 30 && iCSVWgtLF>0. ){
	    wgt_lf_jet30_CSVM += iCSVWgtLF;
	    nowgt_lf_jet30_CSVM += 1.0;
	  }
	  if( csv < csvWPT && pt > 30 && iCSVWgtLF>0. ){
	    wgt_lf_jet30_CSVT += iCSVWgtLF;
	    nowgt_lf_jet30_CSVT += 1.0;
	  }

	  
	  if( csv < csvWPL && pt > 20 && iCSVWgtLF>0. ){
	    wgt_lf_jet20_CSVL += iCSVWgtLF;
	    nowgt_lf_jet20_CSVL += 1.0;
	  }
	  if( csv < csvWPM && pt > 20 && iCSVWgtLF>0. ){
	    wgt_lf_jet20_CSVM += iCSVWgtLF;
	    nowgt_lf_jet20_CSVM += 1.0;
	  }
	  if( csv < csvWPT && pt > 20 && iCSVWgtLF>0. ){
	    wgt_lf_jet20_CSVT += iCSVWgtLF;
	    nowgt_lf_jet20_CSVT += 1.0;
	  }



	  if( cmva < cMVAWPL && pt > 30 && iCMVAWgtLF>0. ){
	    wgt_lf_jet30_CMVAL += iCMVAWgtLF;
	    nowgt_lf_jet30_CMVAL += 1.0;
	  }
	  if( cmva < cMVAWPM && pt > 30 && iCMVAWgtLF>0. ){
	    wgt_lf_jet30_CMVAM += iCMVAWgtLF;
	    nowgt_lf_jet30_CMVAM += 1.0;
	  }
	  if( cmva < cMVAWPT && pt > 30 && iCMVAWgtLF>0. ){
	    wgt_lf_jet30_CMVAT += iCMVAWgtLF;
	    nowgt_lf_jet30_CMVAT += 1.0;
	  }

	  
	  if( cmva < cMVAWPL && pt > 20 && iCMVAWgtLF>0. ){
	    wgt_lf_jet20_CMVAL += iCMVAWgtLF;
	    nowgt_lf_jet20_CMVAL += 1.0;
	  }
	  if( cmva < cMVAWPM && pt > 20 && iCMVAWgtLF>0. ){
	    wgt_lf_jet20_CMVAM += iCMVAWgtLF;
	    nowgt_lf_jet20_CMVAM += 1.0;
	  }
	  if( cmva < cMVAWPT && pt > 20 && iCMVAWgtLF>0. ){
	    wgt_lf_jet20_CMVAT += iCMVAWgtLF;
	    nowgt_lf_jet20_CMVAT += 1.0;
	  }


	}
      }
    if( verbose_ ) std::cout << " ===> test 0.7 " << std::endl;

      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30_CSVL);
      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30_CSVM);
      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30_CSVT);

      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30_CSVL);
      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30_CSVM);
      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30_CSVT);

      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet20_CSVL);
      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet20_CSVM);
      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet20_CSVT);

      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet20_CSVL);
      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet20_CSVM);
      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet20_CSVT);


      h_numEvents_jet30to50_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30to50_CSVL);
      h_numEvents_jet30to50_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30to50_CSVM);
      h_numEvents_jet30to50_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30to50_CSVT);

      h_numEvents_jet30to50_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30to50_CSVL);
      h_numEvents_jet30to50_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30to50_CSVM);
      h_numEvents_jet30to50_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30to50_CSVT);

      
      h_numEvents_jet50to70_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet50to70_CSVL);
      h_numEvents_jet50to70_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet50to70_CSVM);
      h_numEvents_jet50to70_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet50to70_CSVT);

      h_numEvents_jet50to70_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet50to70_CSVL);
      h_numEvents_jet50to70_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet50to70_CSVM);
      h_numEvents_jet50to70_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet50to70_CSVT);

      h_numEvents_jet70to100_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet70to100_CSVL);
      h_numEvents_jet70to100_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet70to100_CSVM);
      h_numEvents_jet70to100_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet70to100_CSVT);

      h_numEvents_jet70to100_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet70to100_CSVL);
      h_numEvents_jet70to100_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet70to100_CSVM);
      h_numEvents_jet70to100_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet70to100_CSVT);

      
      h_numEvents_jet100toInf_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet100toInf_CSVL);
      h_numEvents_jet100toInf_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet100toInf_CSVM);
      h_numEvents_jet100toInf_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet100toInf_CSVT);

      h_numEvents_jet100toInf_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet100toInf_CSVL);
      h_numEvents_jet100toInf_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet100toInf_CSVM);
      h_numEvents_jet100toInf_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet100toInf_CSVT);

      ////

      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30_CMVAL);
      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30_CMVAM);
      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30_CMVAT);

      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30_CMVAL);
      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30_CMVAM);
      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30_CMVAT);

      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet20_CMVAL);
      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet20_CMVAM);
      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet20_CMVAT);

      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet20_CMVAL);
      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet20_CMVAM);
      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet20_CMVAT);


      
      h_numEvents_jet30to50_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30to50_CMVAL);
      h_numEvents_jet30to50_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30to50_CMVAM);
      h_numEvents_jet30to50_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30to50_CMVAT);

      h_numEvents_jet30to50_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30to50_CMVAL);
      h_numEvents_jet30to50_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30to50_CMVAM);
      h_numEvents_jet30to50_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30to50_CMVAT);

      
      h_numEvents_jet50to70_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet50to70_CMVAL);
      h_numEvents_jet50to70_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet50to70_CMVAM);
      h_numEvents_jet50to70_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet50to70_CMVAT);

      h_numEvents_jet50to70_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet50to70_CMVAL);
      h_numEvents_jet50to70_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet50to70_CMVAM);
      h_numEvents_jet50to70_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet50to70_CMVAT);

      h_numEvents_jet70to100_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet70to100_CMVAL);
      h_numEvents_jet70to100_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet70to100_CMVAM);
      h_numEvents_jet70to100_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet70to100_CMVAT);

      h_numEvents_jet70to100_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet70to100_CMVAL);
      h_numEvents_jet70to100_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet70to100_CMVAM);
      h_numEvents_jet70to100_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet70to100_CMVAT);

      
      h_numEvents_jet100toInf_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet100toInf_CMVAL);
      h_numEvents_jet100toInf_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet100toInf_CMVAM);
      h_numEvents_jet100toInf_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet100toInf_CMVAT);

      h_numEvents_jet100toInf_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet100toInf_CMVAL);
      h_numEvents_jet100toInf_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet100toInf_CMVAM);
      h_numEvents_jet100toInf_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet100toInf_CMVAT);

      //mistag
      h_numEvents_lf_jet30_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_lf_jet30_CSVL);
      h_numEvents_lf_jet30_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_lf_jet30_CSVM);
      h_numEvents_lf_jet30_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_lf_jet30_CSVT);

      h_numEvents_lf_jet30_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_lf_jet30_CSVL);
      h_numEvents_lf_jet30_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_lf_jet30_CSVM);
      h_numEvents_lf_jet30_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_lf_jet30_CSVT);

      h_numEvents_lf_jet20_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_lf_jet20_CSVL);
      h_numEvents_lf_jet20_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_lf_jet20_CSVM);
      h_numEvents_lf_jet20_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_lf_jet20_CSVT);

      h_numEvents_lf_jet20_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_lf_jet20_CSVL);
      h_numEvents_lf_jet20_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_lf_jet20_CSVM);
      h_numEvents_lf_jet20_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_lf_jet20_CSVT);



      h_numEvents_lf_jet30_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_lf_jet30_CMVAL);
      h_numEvents_lf_jet30_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_lf_jet30_CMVAM);
      h_numEvents_lf_jet30_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_lf_jet30_CMVAT);

      h_numEvents_lf_jet30_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_lf_jet30_CMVAL);
      h_numEvents_lf_jet30_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_lf_jet30_CMVAM);
      h_numEvents_lf_jet30_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_lf_jet30_CMVAT);

      h_numEvents_lf_jet20_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_lf_jet20_CMVAL);
      h_numEvents_lf_jet20_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_lf_jet20_CMVAM);
      h_numEvents_lf_jet20_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_lf_jet20_CMVAT);

      h_numEvents_lf_jet20_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_lf_jet20_CMVAL);
      h_numEvents_lf_jet20_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_lf_jet20_CMVAM);
      h_numEvents_lf_jet20_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_lf_jet20_CMVAT);


      //// END NEW
    }
    if( verbose_ ) std::cout << " ===> test 2 " << std::endl;

    if( insample<0 ){
      wgt_csv_noSys = 1.0;
      wgt_csv_hf_noSys = 1.0;
      wgt_csv_lf_noSys = 1.0;
      wgt_csv_cf_noSys = 1.0;

      // wgt_csv_noSys_partonFlavourNoPU = 1.0;
      // wgt_csv_hf_noSys_partonFlavourNoPU = 1.0;
      // wgt_csv_lf_noSys_partonFlavourNoPU = 1.0;
      // wgt_csv_cf_noSys_partonFlavourNoPU = 1.0;
    }
    int run = eve->run_;

    //// --------- various weights: PU, topPt, triggerSF, leptonSF...
    // double  wgt_topPtSF = eve->wgt_topPt_;
    double Xsec = mySample_xSec_;//eve->wgt_xs_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    if( maxNentries==1000000 && insample==2300 ) nGen = 670121;
    double lumi = ( intLumi > 0 ) ? intLumi : 10000 ;

    double wgt_lumi = ( insample > 0 ) ? lumi * (Xsec/nGen) : 1;//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights
    wgt_lumi = 1;

    // std::cout << " ************************************************************* " << std::endl;

    // printf("  wgt_CSVL = %.3f, nowgt_CSVL = %.3f, ratioL = %.3f \n", wgt_jet_CSVL, nowgt_jet_CSVL, wgt_jet_CSVL/nowgt_jet_CSVL );
    // printf("  wgt_CSVM = %.3f, nowgt_CSVM = %.3f, ratioM = %.3f \n", wgt_jet_CSVM, nowgt_jet_CSVM, wgt_jet_CSVM/nowgt_jet_CSVM );
    // printf("  wgt_CSVT = %.3f, nowgt_CSVT = %.3f, ratioT = %.3f \n", wgt_jet_CSVT, nowgt_jet_CSVT, wgt_jet_CSVT/nowgt_jet_CSVT );


      
    double wgt = wgt_gen * wgt_lumi * wgt_pu;
    double wgtCSV = wgt * wgt_csv_noSys;

    numEvents_all += 1.;
    numEvents_wgt_gen += wgt_gen;
    numEvents_wgt_lumi += wgt_lumi;
    numEvents_wgt_gen_lumi += wgt_gen * wgt_lumi;
    numEvents_wgt_gen_lumi_pu += wgt_gen * wgt_lumi * wgt_pu;
    numEvents_wgt_gen_lumi_pu_csv += wgt_gen * wgt_lumi * wgt_pu * wgt_csv_noSys;


    //    h_numTruePVs->Fill(numTruePVs,wgt_gen);

    //    h_lheHT_wgt->Fill(lheHT,wgt_gen * wgt_lumi);


    // Skip HT range, if necessary
    // if( useHTbins_ && (insample==2300 || insample==2305 || insample==2400 || insample==2405) ){
    //   if( lheHT > 100 ) continue;
    // }

    //    h_lheHT_wgt_afterHT->Fill(lheHT,wgt_gen * wgt_lumi);

    if( verbose_ ) std::cout << " ===> test 3 " << std::endl;


    ///////////////////
    ////// selections
    ///////////////////


    // require first PV is good
    bool goodFirstVertex = ( eve->goodFirstVertex_ );
    if( !goodFirstVertex ) continue;
    

    bool pass_trigger_ele = false;
    // if( insample<0 ) pass_trigger_ele = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    // else             pass_trigger_ele = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );
    pass_trigger_ele = 1;// || ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    
    bool pass_trigger_mu  = false;
    // if( insample<0 ) pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );
    // else             pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );
    pass_trigger_mu = 1;// || ( eve->pass_HLT_IsoMu20_v_==1 );


    if( insample<0 ){
      pass_trigger_ele = ( insample==-11 && pass_trigger_ele );
      pass_trigger_mu  = ( insample==-13 && pass_trigger_mu  );
    }


    // Pass Trigger selection
    if( !(pass_trigger_mu || pass_trigger_ele) ) continue;

    // require at least one PV
    if( numPVs<1 ) continue;


    // h_numPVs->Fill(numPVs,wgt);


    /////

    int numEle = int( ind_ele.size() );
    int numMu  = int( ind_mu.size() );

    int numEle_loose = int( ind_ele_loose.size() );
    int numMu_loose  = int( ind_mu_loose.size() );

    bool oneEle = ( numEle==1 && numEle_loose==1 && numMu_loose==0 );
    bool oneMu  = ( numMu==1 && numMu_loose==1 && numEle_loose==0 );

    oneEle = ( oneEle && pass_trigger_ele );
    oneMu  = ( oneMu  && pass_trigger_mu );

    bool oneLep = ( oneEle || oneMu );


    bool TwoEle = ( numEle>=1 && numEle_loose==2 && numMu_loose==0 );
    bool TwoMu  = ( numMu>=1 && numMu_loose==2 && numEle_loose==0 );

    TwoEle = ( TwoEle && pass_trigger_ele );
    TwoMu  = ( TwoMu  && pass_trigger_mu );

    bool TwoLep = ( TwoEle || TwoMu );

    bool EleMu = ( (numEle==1 || numMu==1) && numMu_loose==1 && numEle_loose==1 );

    EleMu = ( EleMu && ((numEle==1 && pass_trigger_ele) || (numMu==1 && pass_trigger_mu)) );
    if( verbose_ ) std::cout << " ===> test 3.0.1 " << std::endl;

    // Require at least two leptons

    /// EleMu control region

    // Require at least two leptons
  
    /// end EleMu control region
    if( verbose_ ) std::cout << " ===> test 3.0.4 " << std::endl;

      
    // Require exactly one lepton
    if( !oneLep ) continue;
    if( verbose_ ) std::cout << " ===> test 3.1 " << std::endl;



    int lepInd = ( numEle>0 ) ? ind_ele[0] : ind_mu[0];

    TLorentzVector myLep;
    myLep.SetPtEtaPhiE( eve->lepton_pt_[lepInd], eve->lepton_eta_[lepInd], eve->lepton_phi_[lepInd], eve->lepton_energy_[lepInd] );


    /////
    /////  Loop over systematics
    /////

    for( int iSys=0; iSys<NumSysCat; iSys++ ){
      if( verbose_ ) std::cout << " ===> test 3.3.17 " << std::endl;

      TString sys_name = sys_cat_labels[iSys];

      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_phi = eve->jet_phi_;
      vdouble use_jet_energy = eve->jet_energy_;

      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_cmva = eve->jet_cmva_;
      // vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;




      if( sys_name.Contains("JESUp") ){
	use_jet_pt = eve->jet_JESup_pt_;
	use_jet_eta = eve->jet_JESup_eta_;
	use_jet_phi = eve->jet_JESup_phi_;
	use_jet_energy = eve->jet_JESup_energy_;

	use_jet_csv = eve->jet_JESup_csv_;
	use_jet_cmva = eve->jet_JESup_cmva_;
	// use_jet_pileupJetId_fullDiscriminant = eve->jet_JESup_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;


      }
      else if( sys_name.Contains("JESDown") ){
	use_jet_pt = eve->jet_JESdown_pt_;
	use_jet_eta = eve->jet_JESdown_eta_;
	use_jet_phi = eve->jet_JESdown_phi_;
	use_jet_energy = eve->jet_JESdown_energy_;

	use_jet_csv = eve->jet_JESdown_csv_;
	use_jet_cmva = eve->jet_JESdown_cmva_;
	// use_jet_pileupJetId_fullDiscriminant = eve->jet_JESdown_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;

      }

     if( verbose_ ) std::cout << " ===> test 3.3.18 " << std::endl;

      // Set lepton scale factor equal to 1 for now
      double wgt_lepIdSF = 1.;
      if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = 1.;
      else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = 1.;


      
      // if( insample>=0 ){
      // 	double use_id_pt = std::min(myLep.Pt(), 119.8);
      // 	double use_id_eta = myLep.Eta();
      // 	double use_id_abseta = fabs(use_id_eta);

      // 	double id_sf = 1., iso_sf = 1., hlt_sf = 1.;
	// if( oneEle ){
	//   double use_id_scEta = eve->lepton_eta_[lepInd];

	//   id_sf = h_ele_id_sf->GetBinContent(h_ele_id_sf->FindBin(use_id_scEta,use_id_pt));
	//   iso_sf = h_ele_iso_sf->GetBinContent(h_ele_iso_sf->FindBin(use_id_eta,use_id_pt));
	//   hlt_sf = h_ele_trig_sf->GetBinContent(h_ele_trig_sf->FindBin(use_id_pt,use_id_scEta));
	// }
	// else if( oneMu ){
	//   id_sf = h_mu_id_sf->GetBinContent(h_mu_id_sf->FindBin(use_id_abseta,use_id_pt));
	//   iso_sf = h_mu_iso_sf->GetBinContent(h_mu_iso_sf->FindBin(use_id_abseta,use_id_pt));
	//   hlt_sf = h_mu_trig_sf->GetBinContent(h_mu_trig_sf->FindBin(use_id_abseta,use_id_pt));
	// }

	// wgt_lepIdSF = id_sf * iso_sf * hlt_sf;
	// double hlt_err = ( oneEle ) ? 0.02 : 0.01;
	// if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = id_sf*(1 + 0.01) * iso_sf*(1 + 0.01) * hlt_sf*(1 + hlt_err);
	// else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = id_sf*(1 - 0.01) * iso_sf*(1 - 0.01) * hlt_sf*(1 - hlt_err);

	// if( iSys==0 ){
	//   if( oneEle ) printf("Electron: pT = %.1f, eta = %+.2f, id_sf = %.3f, iso_sf = %.3f, hlt_sf = %.3f, wgt_lepIdSF = %.3f \n", use_id_pt, use_id_eta, id_sf, iso_sf, hlt_sf, wgt_lepIdSF);
	//   else if( oneMu ) printf("Muon: pT = %.1f, eta = %+.2f, id_sf = %.3f, iso_sf = %.3f, hlt_sf = %.3f, wgt_lepIdSF = %.3f \n", use_id_pt, use_id_eta, id_sf, iso_sf, hlt_sf, wgt_lepIdSF);
	// }
      //      }
 
      // PU systematic
      double use_wgt_pu = wgt_pu;
      if( sys_name.Contains("PUUp") )        use_wgt_pu = wgt_pu_up;
      else if( sys_name.Contains("PUDown") ) use_wgt_pu = wgt_pu_down;


     if( verbose_ ) std::cout << " ===> test 3.3.19 " << std::endl;

      // Factorization and normalization scale uncertainty
      double wgt_LHEscale = 1.;
      // if( insample>-1 ){
      // 	double originalXWGTUP = eve->originalXWGTUP_;
      // 	vdouble lhe_event_weights = eve->LHEEvent_weights_;

      // 	if( lhe_event_weights.size()>76 && originalXWGTUP!=0 ){
      // 	  int use_weight_index = 0;

      // 	  if( sys_name.Contains("muFUp") )        use_weight_index = 1;
      // 	  else if( sys_name.Contains("muFDown") ) use_weight_index = 2;
      // 	  else if( sys_name.Contains("muRUp") )   use_weight_index = 3;
      // 	  else if( sys_name.Contains("muRDown") ) use_weight_index = 6;
      // 	  else if( sys_name.Contains("muRmuFUp") )   use_weight_index = 4;
      // 	  else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
      // 	  else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
      // 	  else if( sys_name.Contains("NNPDFUp") ) use_weight_index = 75;
      // 	  else if( sys_name.Contains("NNPDFDown") ) use_weight_index = 13;

      // 	  wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
      // 	}
      // }
     if( verbose_ ) std::cout << " ===> test 3.3.20 " << std::endl;


      //wgt = wgt_gen * wgt_lumi * wgt_pu * wgt_lepIdSF * wgt_LHEscale;
      wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale;


      vdouble jetPts;
      vdouble jetEtas;
      vdouble jetCSVs;
      vdouble jetCMVAs;
      vint    jetFlavors;

      vdouble temp_test_jetPts;
      vdouble temp_test_jetEtas;
      vdouble temp_test_jetCSVs;
      vint    temp_test_jetFlavors;

      vint ind_jet;

      int numJet = 0;
      int numBtag = 0;
      int numBtag_CMVA = 0;

      double HT30=0;
      TLorentzVector sumJet;
      bool firstJet = false;

      // double wgt_csv2=1.;
      // double wgt_csv2_hf=1., wgt_csv2_lf=1., wgt_csv2_cf=1.;

      // double wgt_cmva2=1.;
	    
      if( verbose_ ) std::cout << " ===> test 3.3.21 " << std::endl;

      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	TLorentzVector myJet;
	myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	if( !firstJet ) sumJet = myJet;
	else            sumJet += myJet;

	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];

	if( !(pt>30. && fabs(eta)<2.5) ) continue;

	double dR = myJet.DeltaR(myLep);
	// if( iSys==0 ) h_deltaR_jet_lep->Fill(dR,wgt);

	if( dR < 0.4 ) continue;

	double csv = use_jet_csv[iJet];
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	double cmva = use_jet_cmva[iJet];

	int flavor = use_jet_hadronFlavour[iJet];

	ind_jet.push_back(iJet);

	HT30 += pt;

	temp_test_jetPts.clear();
	temp_test_jetEtas.clear();
	temp_test_jetCSVs.clear();
	temp_test_jetFlavors.clear();

	temp_test_jetPts.push_back(pt);
	temp_test_jetEtas.push_back(eta);
	temp_test_jetCSVs.push_back(csv);
	temp_test_jetFlavors.push_back(flavor);


	jetPts.push_back(pt);
	jetEtas.push_back(eta);
	jetCSVs.push_back(csv);
	jetCMVAs.push_back(cmva);
	jetFlavors.push_back(flavor);

	numJet++;
	if( csv > csvWPM ){
	  numBtag++;
	}
	if( cmva > cMVAWPM ){
	  numBtag_CMVA++;
	}

	if( pt > 1000 ) pt = 999.;



	/////	
	double temp_test_wgt_csv_hf=1, temp_test_wgt_csv_lf=1, temp_test_wgt_csv_cf=1;
	double temp_test_wgt_csv = get_csv_wgt(temp_test_jetPts, temp_test_jetEtas, temp_test_jetCSVs, temp_test_jetFlavors, iSys, 
					       temp_test_wgt_csv_hf, temp_test_wgt_csv_lf, temp_test_wgt_csv_cf);

	if( iSys==0 ){
	  if( abs(flavor)==5 ){
	    h_b_csv_wgt->Fill(csv,temp_test_wgt_csv);
	  }
	  else if( abs(flavor)==4 ){
	    h_c_csv_wgt->Fill(csv,temp_test_wgt_csv);
	  }
	  else {
	    h_l_csv_wgt->Fill(csv,temp_test_wgt_csv);
	  }
	}

	// if( sys_name.Contains("JESUp") && false ){
	//   printf("  iJet = %d: pt = %4.1f, eta = %+.2f, flavor = %+d, CSVv2 = %.2f: weight csv file = %.3f, weight hist file = %.3f \n",
	// 	 iJet, pt, eta, flavor, csv, my_jet_sf, temp_test_wgt_csv);
	// }


      }


      // Calculate CSV weight
      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
      double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      wgtCSV = wgt * wgt_csv;
    if( verbose_ ) std::cout << " ===> test 3.3.1 " << std::endl;


      // Calculate CMVA weight
      double wgt_cmva_hf, wgt_cmva_lf, wgt_cmva_cf;
      double wgt_cmva = ( insample<0 ) ? 1 : get_cmva_wgt(jetPts, jetEtas, jetCMVAs, jetFlavors, iSys, wgt_cmva_hf, wgt_cmva_lf, wgt_cmva_cf);
    if( verbose_ ) std::cout << " ===> test 3.3.4 " << std::endl;

      double wgtCMVA = wgt * wgt_cmva;

      // if( insample<0 ) wgt_cmva2 = 1;
      // double wgtCMVA2 = wgt * wgt_cmva2;

      
      //if( iSys==0 ) printf("\t total event weight: csv file = %.3f, hist file = %.3f, btv csv file = %.3f \n", wgt_csv2, wgt_csv, wgt_csv3);

      //std::cout << "\t main: iSys = " << iSys << ",\t get_csv_wgt wgt_csv_lf = " << wgt_csv_lf << std::endl;


      h_hf_wgtCSV_perSys->Fill(iSys,wgt_csv_hf);
      h_lf_wgtCSV_perSys->Fill(iSys,wgt_csv_lf);
      h_cf_wgtCSV_perSys->Fill(iSys,wgt_csv_cf);



      h_numEvents_perSys->Fill(iSys,wgt);
      h_numEvents_perSys_wgtCSV->Fill(iSys,wgtCSV);
      if( iSys==0 ) h_numEvents_perSys_wgtCSV->Fill(-1,wgt);

      // h_numEvents_perSys_wgtCSV2->Fill(iSys,wgtCSV2);
      // if( iSys==0 ) h_numEvents_perSys_wgtCSV2->Fill(-1,wgt);

      h_numEvents_perSys_wgtCMVA->Fill(iSys,wgtCMVA);
      if( iSys==0 ) h_numEvents_perSys_wgtCMVA->Fill(-1,wgt);

      // h_numEvents_perSys_wgtCMVA2->Fill(iSys,wgtCMVA2);
      // if( iSys==0 ) h_numEvents_perSys_wgtCMVA2->Fill(-1,wgt);


      // int fill_numJet  = std::min( MaxNjet-1, numJet );
      // int fill_numBtag = std::min( MaxNbtag-1, numBtag );
      // int fill_numBtag_CMVA = std::min( MaxNbtag-1, numBtag_CMVA );





      // Require at least four jets
      if( !(numJet>=4) ) continue;




    if( verbose_ ) std::cout << " ===> test 3.3.6 " << std::endl;


      // double fill_pfMETNoHF = std::min( metmax-0.001, pfMETNoHF);

      // // Fill MET histograms

    if( verbose_ ) std::cout << " ===> test 3.3.7 " << std::endl;


    if( verbose_ ) std::cout << " ===> test 3.3.9 " << std::endl;

      
      // Require at least two b-tags
      if( !(numBtag>=2) ) continue;



      // h_numEvents_perSys_wgtCSV_4j2t->Fill(iSys,wgtCSV);
      // if( iSys==0 ) h_numEvents_perSys_wgtCSV_4j2t->Fill(-1,wgt);

      // if( numJet>=6 && numBtag>=4 ){
      // 	h_numEvents_perSys_wgtCSV_6j4t->Fill(iSys,wgtCSV);
      // 	if( iSys==0 ) h_numEvents_perSys_wgtCSV_6j4t->Fill(-1,wgt);
      // }

      
      // h_numEvents_perSys_wgtCSV2_4j2t->Fill(iSys,wgtCSV2);
      // if( iSys==0 ) h_numEvents_perSys_wgtCSV2_4j2t->Fill(-1,wgt);

      // if( numJet>=6 && numBtag>=4 ){
      // 	h_numEvents_perSys_wgtCSV2_6j4t->Fill(iSys,wgtCSV2);
      // 	if( iSys==0 ) h_numEvents_perSys_wgtCSV2_6j4t->Fill(-1,wgt);
      // }
      
      numEvents_wgt_gen_lumi_pu_4j2t += wgt_gen * wgt_lumi * wgt_pu;
      numEvents_wgt_gen_lumi_pu_csv_4j2t += wgt_gen * wgt_lumi * wgt_pu * wgt_csv_noSys;




      ////////////

    if( verbose_ ) std::cout << " ===> test 3.3.10 " << std::endl;



      // std::sort(jetPts.begin(), jetPts.end());
      // std::reverse(jetPts.begin(), jetPts.end());

      // if( jetPts.size()>=1 ) h_jet_1_pt->Fill(jetPts[0],wgt);
      // if( jetPts.size()>=2 ) h_jet_2_pt->Fill(jetPts[1],wgt);
      // if( jetPts.size()>=3 ) h_jet_3_pt->Fill(jetPts[2],wgt);
      // if( jetPts.size()>=4 ) h_jet_4_pt->Fill(jetPts[3],wgt);



    if( verbose_ ) std::cout << " ===> test 3.3.11 " << std::endl;


      //
      // Fill MET
      //



      //MHT


    if( verbose_ ) std::cout << " ===> test 3.3.12 " << std::endl;

    } // End loop over systematics

    if( verbose_ ) std::cout << " ===> test 3.3.13 " << std::endl;

  } // end loop over events

  std::cout << "**************************************************************" << std::endl;
  std::cout << "\t Number of all events = " << numEvents_all << std::endl;
  std::cout << "\t Number of gen weighted events = " << numEvents_wgt_gen << std::endl;
  std::cout << "\t Number of lumi weighted events = " << numEvents_wgt_lumi << std::endl;
  std::cout << "\t Number of gen * lumi weighted events = " << numEvents_wgt_gen_lumi << std::endl;
  std::cout << "\t Number of gen * lumi * PU weighted events = " << numEvents_wgt_gen_lumi_pu << std::endl;
  std::cout << "\t Number of gen * lumi * PU * CSV weighted events = " << numEvents_wgt_gen_lumi_pu_csv << std::endl;
  std::cout << "\t Number of gen * lumi * PU weighted events (>=4j,>=2b) = " << numEvents_wgt_gen_lumi_pu_4j2t << std::endl;
  std::cout << "\t Number of gen * lumi * PU * CSV weighted events (>=4j,>=2b) = " << numEvents_wgt_gen_lumi_pu_csv_4j2t  << std::endl;
  std::cout << "**************************************************************" << std::endl;



  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();

}


double reweightPU( int nPU, int iSys ){

  double PUscale[50];

  if( nPU>49 ) nPU = 49;

  if( iSys==0 ){
  PUscale[0] = 0.544141;
  PUscale[1] = 0.683152;
  PUscale[2] = 1.0588;
  PUscale[3] = 1.36648;
  PUscale[4] = 1.62214;
  PUscale[5] = 1.93965;
  PUscale[6] = 1.45916;
  PUscale[7] = 1.29301;
  PUscale[8] = 1.39344;
  PUscale[9] = 1.37393;
  PUscale[10] = 1.26638;
  PUscale[11] = 1.1637;
  PUscale[12] = 1.05777;
  PUscale[13] = 0.901478;
  PUscale[14] = 0.695434;
  PUscale[15] = 0.48665;
  PUscale[16] = 0.323436;
  PUscale[17] = 0.233179;
  PUscale[18] = 0.204042;
  PUscale[19] = 0.179695;
  PUscale[20] = 0.117906;
  PUscale[21] = 0.0535705;
  PUscale[22] = 0.0179761;
  PUscale[23] = 0.00486845;
  PUscale[24] = 0.00120473;
  PUscale[25] = 0.000306356;
  PUscale[26] = 8.44097e-05;
  PUscale[27] = 2.35712e-05;
  PUscale[28] = 6.33309e-06;
  PUscale[29] = 1.62779e-06;
  PUscale[30] = 3.92723e-07;
  PUscale[31] = 8.94453e-08;
  PUscale[32] = 1.90178e-08;
  PUscale[33] = 3.84203e-09;
  PUscale[34] = 7.38767e-10;
  PUscale[35] = 1.28456e-10;
  PUscale[36] = 2.08504e-11;
  PUscale[37] = 3.31182e-12;
  PUscale[38] = 4.63316e-13;
  PUscale[39] = 4.72854e-14;
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
  }
  else if( iSys==1 ){
  PUscale[0] = 0.458261;
  PUscale[1] = 0.626711;
  PUscale[2] = 0.971785;
  PUscale[3] = 1.21776;
  PUscale[4] = 1.34764;
  PUscale[5] = 1.4248;
  PUscale[6] = 0.922526;
  PUscale[7] = 0.744979;
  PUscale[8] = 0.913723;
  PUscale[9] = 1.06132;
  PUscale[10] = 1.11114;
  PUscale[11] = 1.12292;
  PUscale[12] = 1.13425;
  PUscale[13] = 1.10829;
  PUscale[14] = 1.00558;
  PUscale[15] = 0.83577;
  PUscale[16] = 0.639871;
  PUscale[17] = 0.469319;
  PUscale[18] = 0.361301;
  PUscale[19] = 0.306711;
  PUscale[20] = 0.237734;
  PUscale[21] = 0.142901;
  PUscale[22] = 0.0656738;
  PUscale[23] = 0.0238299;
  PUscale[24] = 0.00722283;
  PUscale[25] = 0.00197294;
  PUscale[26] = 0.000550982;
  PUscale[27] = 0.000165108;
  PUscale[28] = 5.14307e-05;
  PUscale[29] = 1.60005e-05;
  PUscale[30] = 4.76288e-06;
  PUscale[31] = 1.35321e-06;
  PUscale[32] = 3.62213e-07;
  PUscale[33] = 9.29401e-08;
  PUscale[34] = 2.28987e-08;
  PUscale[35] = 5.14682e-09;
  PUscale[36] = 1.08953e-09;
  PUscale[37] = 2.27672e-10;
  PUscale[38] = 4.24244e-11;
  PUscale[39] = 6.6564e-12;
  PUscale[40] = 1.21726e-12;
  PUscale[41] = 1.25024e-13;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }
  else if( iSys==-1 ){
  PUscale[0] = 0.637352;
  PUscale[1] = 0.745524;
  PUscale[2] = 1.16977;
  PUscale[3] = 1.55999;
  PUscale[4] = 2.01896;
  PUscale[5] = 2.79883;
  PUscale[6] = 2.50005;
  PUscale[7] = 2.22714;
  PUscale[8] = 2.00474;
  PUscale[9] = 1.67907;
  PUscale[10] = 1.37446;
  PUscale[11] = 1.13578;
  PUscale[12] = 0.899819;
  PUscale[13] = 0.646461;
  PUscale[14] = 0.412783;
  PUscale[15] = 0.243661;
  PUscale[16] = 0.153684;
  PUscale[17] = 0.125858;
  PUscale[18] = 0.116543;
  PUscale[19] = 0.0860169;
  PUscale[20] = 0.0414331;
  PUscale[21] = 0.0133115;
  PUscale[22] = 0.00326222;
  PUscale[23] = 0.000723524;
  PUscale[24] = 0.000169442;
  PUscale[25] = 4.20509e-05;
  PUscale[26] = 1.04015e-05;
  PUscale[27] = 2.41401e-06;
  PUscale[28] = 5.19674e-07;
  PUscale[29] = 1.0525e-07;
  PUscale[30] = 1.97832e-08;
  PUscale[31] = 3.47419e-09;
  PUscale[32] = 5.63765e-10;
  PUscale[33] = 8.60403e-11;
  PUscale[34] = 1.23701e-11;
  PUscale[35] = 1.59276e-12;
  PUscale[36] = 1.89108e-13;
  PUscale[37] = 1.98468e-14;
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
  }

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
    for( int iPt=0; iPt<4; iPt++ ){
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





void fillCMVAhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_cmva_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_cmva_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_cmva_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_cmva_suffix_hf = "final";
    TString syst_cmva_suffix_c = "final";
    TString syst_cmva_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_cmva_suffix_hf = "final_JESUp"; syst_cmva_suffix_lf = "final_JESUp";
      syst_cmva_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_cmva_suffix_hf = "final_JESDown"; syst_cmva_suffix_lf = "final_JESDown";
      syst_cmva_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_cmva_suffix_hf = "final_LFUp"; syst_cmva_suffix_lf = "final_HFUp";
      syst_cmva_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_cmva_suffix_hf = "final_LFDown"; syst_cmva_suffix_lf = "final_HFDown";
      syst_cmva_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_cmva_suffix_hf = "final_Stats1Up"; syst_cmva_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_cmva_suffix_hf = "final_Stats1Down"; syst_cmva_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_cmva_suffix_hf = "final_Stats2Up"; syst_cmva_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_cmva_suffix_hf = "final_Stats2Down"; syst_cmva_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_cmva_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_cmva_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_cmva_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_cmva_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_cmva_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_cmva_suffix_lf.Data()) );
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
    else if (pt >=100)          iPt = 3;
    //else if (pt >=100 && pt<160) iPt = 4;
    //else if (pt >=160 && pt<10000) iPt = 5;

    ////
    if(abs(flavor) == 5 || abs(flavor) == 4){
      if (pt >=19.99 && pt<30) iPt = 0;
      else if (pt >=30 && pt<50) iPt = 1;
      else if (pt >=50 && pt<70) iPt = 2;
      else if (pt >=70 && pt<100) iPt = 3;
      else if (pt >=100 && pt<160) iPt = 4;
      else if (pt >=160 && pt<10000) iPt = 4; //4
    }

    if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.51 ) iEta = 2;

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
      // if( iSysLF==0 ) printf(" iSysLF = %d, iJet = %d,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     iSysLF, iJet, flavor, pt, jetAbsEta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  //std::cout << "\t iSys = " << iSys << " total csvWgtLF = " << csvWgtLF << std::endl;

  return csvWgtTotal;
}





double get_cmva_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCMVA, vint jetFlavor, int iSys, double &cmvaWgtHF, double &cmvaWgtLF, double &cmvaWgtCF ){

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

  double cmvaWgthf = 1.;
  double cmvaWgtC  = 1.;
  double cmvaWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPt.size()); iJet++ ){

    double cmva = jetCMVA[iJet];
    double pt = jetPt[iJet];
    double jetAbsEta = fabs(jetEta[iJet]);
    int flavor = jetFlavor[iJet];

    int iPt = -1; int iEta = -1;
    if (pt >=19.99 && pt<30) iPt = 0;
    else if (pt >=30 && pt<40) iPt = 1;
    else if (pt >=40 && pt<60) iPt = 2;
    else if (pt >=60 && pt<100) iPt = 3;
    else if (pt >=100)          iPt = 3;
    //else if (pt >=100 && pt<160) iPt = 4;
    //else if (pt >=160 && pt<10000) iPt = 5;

    ////
    if(abs(flavor) == 5 || abs(flavor) == 4){
      if (pt >=19.99 && pt<30) iPt = 0;
      else if (pt >=30 && pt<50) iPt = 1;
      else if (pt >=50 && pt<70) iPt = 2;
      else if (pt >=70 && pt<100) iPt = 3;
      else if (pt >=100 && pt<160) iPt = 4;
      else if (pt >=160 && pt<10000) iPt = 4; //4
    }

    if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.51 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;

    //std::cout << " flavor = " << flavor << ", cmva = " << cmva << ", iPt = " << iPt << ", iEta = " << iEta << ", iSysHF = " << iSysHF << ", iSysC = " << iSysC << ", iSysLF = " << iSysLF << std::endl;

    if (abs(flavor) == 5 ){
      int useCMVABin = (cmva>=0.) ? h_cmva_wgt_hf[iSysHF][iPt]->FindBin(cmva) : 1;
      double iCMVAWgtHF = h_cmva_wgt_hf[iSysHF][iPt]->GetBinContent(useCMVABin);
      if( iCMVAWgtHF!=0 ) cmvaWgthf *= iCMVAWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, cmva, iCMVAWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCMVABin = (cmva>=0.) ? c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva) : 1;
      double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);
      if( iCMVAWgtC!=0 ) cmvaWgtC *= iCMVAWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      // 			    flavor, pt, jetAbsEta, cmva, iCMVAWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCMVABin = (cmva>=0.) ? h_cmva_wgt_lf[iSysLF][iPt][iEta]->FindBin(cmva) : 1;
      double iCMVAWgtLF = h_cmva_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCMVABin);
      if( iCMVAWgtLF!=0 ) cmvaWgtlf *= iCMVAWgtLF;
      // if( iSysLF==0 ) printf(" iSysLF = %d, iJet = %d,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      //  			     iSysLF, iJet, flavor, pt, jetAbsEta, cmva, iCMVAWgtLF );
    }
  }

  double cmvaWgtTotal = cmvaWgthf * cmvaWgtC * cmvaWgtlf;

  cmvaWgtHF = cmvaWgthf;
  cmvaWgtLF = cmvaWgtlf;
  cmvaWgtCF = cmvaWgtC;

  //std::cout << "\t iSys = " << iSys << " total cmvaWgtLF = " << cmvaWgtLF << std::endl;

  return cmvaWgtTotal;
}



// void fillCSVEffhistos(TFile* file){

//   h_a_jet_pt_eta_all_ = (TH2D*)file->Get("h_a_jet_pt_eta_all")->Clone("h_a_jet_pt_eta_all_temp");
//   h_b_jet_pt_eta_all_ = (TH2D*)file->Get("h_b_jet_pt_eta_all")->Clone("h_b_jet_pt_eta_all_temp");
//   h_c_jet_pt_eta_all_ = (TH2D*)file->Get("h_c_jet_pt_eta_all")->Clone("h_c_jet_pt_eta_all_temp");
//   h_l_jet_pt_eta_all_ = (TH2D*)file->Get("h_l_jet_pt_eta_all")->Clone("h_l_jet_pt_eta_all_temp");

//   h_a_jet_pt_eta_eff_ = (TH2D*)file->Get("h_a_jet_pt_eta_csvM")->Clone("h_a_jet_pt_eta_eff_temp");
//   h_b_jet_pt_eta_eff_ = (TH2D*)file->Get("h_b_jet_pt_eta_csvM")->Clone("h_b_jet_pt_eta_eff_temp");
//   h_c_jet_pt_eta_eff_ = (TH2D*)file->Get("h_c_jet_pt_eta_csvM")->Clone("h_c_jet_pt_eta_eff_temp");
//   h_l_jet_pt_eta_eff_ = (TH2D*)file->Get("h_l_jet_pt_eta_csvM")->Clone("h_l_jet_pt_eta_eff_temp");

//   h_a_jet_pt_eta_eff_->Divide(h_a_jet_pt_eta_all_);
//   h_b_jet_pt_eta_eff_->Divide(h_b_jet_pt_eta_all_);
//   h_c_jet_pt_eta_eff_->Divide(h_c_jet_pt_eta_all_);
//   h_l_jet_pt_eta_eff_->Divide(h_l_jet_pt_eta_all_);
// }




