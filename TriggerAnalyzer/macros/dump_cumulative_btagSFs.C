//Root includes                                   
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TAxis.h"
#include "TKey.h"
#include "TList.h"

void dump_cumulative_btagSFs(bool isCSV = true, bool isHF = true, TString inputFileName  = "file.root") {

  TFile *histFile = TFile::Open(inputFileName);

  TString btagger = ( isCSV ) ? "CSV" : "CMVA";
  TString btaggerStr = ( isCSV ) ? "DeepCSV" : "CSVv2";
  TString flavor = (isHF) ? "" : "lf_";
  TString flavorStr = (isHF) ? "HF" : "LF";
  
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


  
  
  TH1D* h_numEvents_jet30_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet30to50_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30to50_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet50to70_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet50to70_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet70to100_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet70to100_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet100toInf_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet100toInf_nowgtCMVA_WP[numCSVSys];

  double sf_jet30_csvL[numCSVSys];
  double sf_jet30_csvM[numCSVSys];
  double sf_jet30_csvT[numCSVSys];
  
  double sf_jet20_csvL[numCSVSys];
  double sf_jet20_csvM[numCSVSys];
  double sf_jet20_csvT[numCSVSys];

  double sf_jet30to50_cmvaL[numCSVSys];
  double sf_jet30to50_cmvaM[numCSVSys];
  double sf_jet30to50_cmvaT[numCSVSys];
  
  double sf_jet50to70_cmvaL[numCSVSys];
  double sf_jet50to70_cmvaM[numCSVSys];
  double sf_jet50to70_cmvaT[numCSVSys];
  
  double sf_jet70to100_cmvaL[numCSVSys];
  double sf_jet70to100_cmvaM[numCSVSys];
  double sf_jet70to100_cmvaT[numCSVSys];
  
  double sf_jet100toInf_cmvaL[numCSVSys];
  double sf_jet100toInf_cmvaM[numCSVSys];
  double sf_jet100toInf_cmvaT[numCSVSys];

  double sf_jet30_csvL_err_up2 = 0;
  double sf_jet30_csvM_err_up2 = 0;
  double sf_jet30_csvT_err_up2 = 0;

  double sf_jet30_csvL_err_down2 = 0;
  double sf_jet30_csvM_err_down2 = 0;
  double sf_jet30_csvT_err_down2 = 0;


  double sf_jet20_csvL_err_up2 = 0;
  double sf_jet20_csvM_err_up2 = 0;
  double sf_jet20_csvT_err_up2 = 0;
  
  double sf_jet20_csvL_err_down2 = 0;
  double sf_jet20_csvM_err_down2 = 0;
  double sf_jet20_csvT_err_down2 = 0;




  double sf_jet30to50_cmvaL_err_up2 = 0;
  double sf_jet30to50_cmvaM_err_up2 = 0;
  double sf_jet30to50_cmvaT_err_up2 = 0;

  double sf_jet30to50_cmvaL_err_down2 = 0;
  double sf_jet30to50_cmvaM_err_down2 = 0;
  double sf_jet30to50_cmvaT_err_down2 = 0;


  double sf_jet50to70_cmvaL_err_up2 = 0;
  double sf_jet50to70_cmvaM_err_up2 = 0;
  double sf_jet50to70_cmvaT_err_up2 = 0;

  double sf_jet50to70_cmvaL_err_down2 = 0;
  double sf_jet50to70_cmvaM_err_down2 = 0;
  double sf_jet50to70_cmvaT_err_down2 = 0;


  double sf_jet70to100_cmvaL_err_up2 = 0;
  double sf_jet70to100_cmvaM_err_up2 = 0;
  double sf_jet70to100_cmvaT_err_up2 = 0;

  double sf_jet70to100_cmvaL_err_down2 = 0;
  double sf_jet70to100_cmvaM_err_down2 = 0;
  double sf_jet70to100_cmvaT_err_down2 = 0;


  double sf_jet100toInf_cmvaL_err_up2 = 0;
  double sf_jet100toInf_cmvaM_err_up2 = 0;
  double sf_jet100toInf_cmvaT_err_up2 = 0;

  double sf_jet100toInf_cmvaL_err_down2 = 0;
  double sf_jet100toInf_cmvaM_err_down2 = 0;
  double sf_jet100toInf_cmvaT_err_down2 = 0;

  
  for( int iSys=0; iSys<numCSVSys; iSys++ ){
    int useSys = csv_systematics[iSys];

    h_numEvents_jet30_wgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_%sjet30_wgt%s_WP_iSys_%d",flavor.Data(),btagger.Data(),useSys));
    h_numEvents_jet30_nowgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_%sjet30_nowgt%s_WP_iSys_%d",flavor.Data(),btagger.Data(),useSys));

    h_numEvents_jet20_wgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_%sjet20_wgt%s_WP_iSys_%d",flavor.Data(),btagger.Data(),useSys));
    h_numEvents_jet20_nowgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_%sjet20_nowgt%s_WP_iSys_%d",flavor.Data(),btagger.Data(),useSys));


    h_numEvents_jet30to50_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30to50_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet30to50_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30to50_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));

    h_numEvents_jet50to70_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet50to70_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet50to70_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet50to70_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));

    h_numEvents_jet70to100_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet70to100_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet70to100_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet70to100_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));

    h_numEvents_jet100toInf_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet100toInf_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet100toInf_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet100toInf_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));




    h_numEvents_jet30_wgtCSV_WP[iSys]->Divide(h_numEvents_jet30_nowgtCSV_WP[iSys]);
    h_numEvents_jet20_wgtCSV_WP[iSys]->Divide(h_numEvents_jet20_nowgtCSV_WP[iSys]);

    
    h_numEvents_jet30to50_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet30to50_nowgtCMVA_WP[iSys]);
    h_numEvents_jet50to70_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet50to70_nowgtCMVA_WP[iSys]);
    h_numEvents_jet70to100_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet70to100_nowgtCMVA_WP[iSys]);
    h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet100toInf_nowgtCMVA_WP[iSys]);

    if( iSys==0 ){
      sf_jet30_csvL[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(1);
      sf_jet30_csvM[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(2);
      sf_jet30_csvT[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(3);

      sf_jet20_csvL[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(1);
      sf_jet20_csvM[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(2);
      sf_jet20_csvT[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(3);


      sf_jet30to50_cmvaL[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet30to50_cmvaM[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet30to50_cmvaT[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(3);

      sf_jet50to70_cmvaL[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet50to70_cmvaM[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet50to70_cmvaT[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(3);

      sf_jet70to100_cmvaL[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet70to100_cmvaM[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet70to100_cmvaT[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(3);

      sf_jet100toInf_cmvaL[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet100toInf_cmvaM[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet100toInf_cmvaT[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(3);

    }
    if( iSys>0 ){
      sf_jet30_csvL[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(1) - sf_jet30_csvL[0];
      sf_jet30_csvM[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(2) - sf_jet30_csvM[0];
      sf_jet30_csvT[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(3) - sf_jet30_csvT[0];

      sf_jet20_csvL[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(1) - sf_jet20_csvL[0];
      sf_jet20_csvM[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(2) - sf_jet20_csvM[0];
      sf_jet20_csvT[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(3) - sf_jet20_csvT[0];



      sf_jet30to50_cmvaL[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet30to50_cmvaL[0];
      sf_jet30to50_cmvaM[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet30to50_cmvaM[0];
      sf_jet30to50_cmvaT[iSys] = h_numEvents_jet30to50_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet30to50_cmvaT[0];

      sf_jet50to70_cmvaL[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet50to70_cmvaL[0];
      sf_jet50to70_cmvaM[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet50to70_cmvaM[0];
      sf_jet50to70_cmvaT[iSys] = h_numEvents_jet50to70_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet50to70_cmvaT[0];

      sf_jet70to100_cmvaL[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet70to100_cmvaL[0];
      sf_jet70to100_cmvaM[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet70to100_cmvaM[0];
      sf_jet70to100_cmvaT[iSys] = h_numEvents_jet70to100_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet70to100_cmvaT[0];

      sf_jet100toInf_cmvaL[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet100toInf_cmvaL[0];
      sf_jet100toInf_cmvaM[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet100toInf_cmvaM[0];
      sf_jet100toInf_cmvaT[iSys] = h_numEvents_jet100toInf_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet100toInf_cmvaT[0];

      
      if( sf_jet30_csvL[iSys] > 0 ) sf_jet30_csvL_err_up2   += sf_jet30_csvL[iSys] * sf_jet30_csvL[iSys];
      else                          sf_jet30_csvL_err_down2 += sf_jet30_csvL[iSys] * sf_jet30_csvL[iSys];

      if( sf_jet30_csvM[iSys] > 0 ) sf_jet30_csvM_err_up2   += sf_jet30_csvM[iSys] * sf_jet30_csvM[iSys];
      else                          sf_jet30_csvM_err_down2 += sf_jet30_csvM[iSys] * sf_jet30_csvM[iSys];

      if( sf_jet30_csvT[iSys] > 0 ) sf_jet30_csvT_err_up2   += sf_jet30_csvT[iSys] * sf_jet30_csvT[iSys];
      else                          sf_jet30_csvT_err_down2 += sf_jet30_csvT[iSys] * sf_jet30_csvT[iSys];

      
      if( sf_jet20_csvL[iSys] > 0 ) sf_jet20_csvL_err_up2   += sf_jet20_csvL[iSys] * sf_jet20_csvL[iSys];
      else                          sf_jet20_csvL_err_down2 += sf_jet20_csvL[iSys] * sf_jet20_csvL[iSys];

      if( sf_jet20_csvM[iSys] > 0 ) sf_jet20_csvM_err_up2   += sf_jet20_csvM[iSys] * sf_jet20_csvM[iSys];
      else                          sf_jet20_csvM_err_down2 += sf_jet20_csvM[iSys] * sf_jet20_csvM[iSys];

      if( sf_jet20_csvT[iSys] > 0 ) sf_jet20_csvT_err_up2   += sf_jet20_csvT[iSys] * sf_jet20_csvT[iSys];
      else                          sf_jet20_csvT_err_down2 += sf_jet20_csvT[iSys] * sf_jet20_csvT[iSys];

      

      
      if( sf_jet30to50_cmvaL[iSys] > 0 ) sf_jet30to50_cmvaL_err_up2   += sf_jet30to50_cmvaL[iSys] * sf_jet30to50_cmvaL[iSys];
      else                               sf_jet30to50_cmvaL_err_down2 += sf_jet30to50_cmvaL[iSys] * sf_jet30to50_cmvaL[iSys];

      if( sf_jet30to50_cmvaM[iSys] > 0 ) sf_jet30to50_cmvaM_err_up2   += sf_jet30to50_cmvaM[iSys] * sf_jet30to50_cmvaM[iSys];
      else                               sf_jet30to50_cmvaM_err_down2 += sf_jet30to50_cmvaM[iSys] * sf_jet30to50_cmvaM[iSys];

      if( sf_jet30to50_cmvaT[iSys] > 0 ) sf_jet30to50_cmvaT_err_up2   += sf_jet30to50_cmvaT[iSys] * sf_jet30to50_cmvaT[iSys];
      else                               sf_jet30to50_cmvaT_err_down2 += sf_jet30to50_cmvaT[iSys] * sf_jet30to50_cmvaT[iSys];

      
      if( sf_jet50to70_cmvaL[iSys] > 0 ) sf_jet50to70_cmvaL_err_up2   += sf_jet50to70_cmvaL[iSys] * sf_jet50to70_cmvaL[iSys];
      else                                sf_jet50to70_cmvaL_err_down2 += sf_jet50to70_cmvaL[iSys] * sf_jet50to70_cmvaL[iSys];

      if( sf_jet50to70_cmvaM[iSys] > 0 ) sf_jet50to70_cmvaM_err_up2   += sf_jet50to70_cmvaM[iSys] * sf_jet50to70_cmvaM[iSys];
      else                                sf_jet50to70_cmvaM_err_down2 += sf_jet50to70_cmvaM[iSys] * sf_jet50to70_cmvaM[iSys];

      if( sf_jet50to70_cmvaT[iSys] > 0 ) sf_jet50to70_cmvaT_err_up2   += sf_jet50to70_cmvaT[iSys] * sf_jet50to70_cmvaT[iSys];
      else                                sf_jet50to70_cmvaT_err_down2 += sf_jet50to70_cmvaT[iSys] * sf_jet50to70_cmvaT[iSys];


      if( sf_jet70to100_cmvaL[iSys] > 0 ) sf_jet70to100_cmvaL_err_up2   += sf_jet70to100_cmvaL[iSys] * sf_jet70to100_cmvaL[iSys];
      else                                sf_jet70to100_cmvaL_err_down2 += sf_jet70to100_cmvaL[iSys] * sf_jet70to100_cmvaL[iSys];

      if( sf_jet70to100_cmvaM[iSys] > 0 ) sf_jet70to100_cmvaM_err_up2   += sf_jet70to100_cmvaM[iSys] * sf_jet70to100_cmvaM[iSys];
      else                                sf_jet70to100_cmvaM_err_down2 += sf_jet70to100_cmvaM[iSys] * sf_jet70to100_cmvaM[iSys];

      if( sf_jet70to100_cmvaT[iSys] > 0 ) sf_jet70to100_cmvaT_err_up2   += sf_jet70to100_cmvaT[iSys] * sf_jet70to100_cmvaT[iSys];
      else                                sf_jet70to100_cmvaT_err_down2 += sf_jet70to100_cmvaT[iSys] * sf_jet70to100_cmvaT[iSys];

      
      if( sf_jet100toInf_cmvaL[iSys] > 0 ) sf_jet100toInf_cmvaL_err_up2   += sf_jet100toInf_cmvaL[iSys] * sf_jet100toInf_cmvaL[iSys];
      else                                 sf_jet100toInf_cmvaL_err_down2 += sf_jet100toInf_cmvaL[iSys] * sf_jet100toInf_cmvaL[iSys];

      if( sf_jet100toInf_cmvaM[iSys] > 0 ) sf_jet100toInf_cmvaM_err_up2   += sf_jet100toInf_cmvaM[iSys] * sf_jet100toInf_cmvaM[iSys];
      else                                 sf_jet100toInf_cmvaM_err_down2 += sf_jet100toInf_cmvaM[iSys] * sf_jet100toInf_cmvaM[iSys];

      if( sf_jet100toInf_cmvaT[iSys] > 0 ) sf_jet100toInf_cmvaT_err_up2   += sf_jet100toInf_cmvaT[iSys] * sf_jet100toInf_cmvaT[iSys];
      else                                 sf_jet100toInf_cmvaT_err_down2 += sf_jet100toInf_cmvaT[iSys] * sf_jet100toInf_cmvaT[iSys];

    }
  }


  
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT > 30 " << std::endl;
  std::cout << " " << std::endl;
  if(isHF){
  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvL[0], sqrt(sf_jet30_csvL_err_up2), sqrt(sf_jet30_csvL_err_down2), sf_jet30_csvL[1], sf_jet30_csvL[2], sf_jet30_csvL[3], sf_jet30_csvL[4],
	 sf_jet30_csvL[7], sf_jet30_csvL[8], sf_jet30_csvL[9], sf_jet30_csvL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvM[0], sqrt(sf_jet30_csvM_err_up2), sqrt(sf_jet30_csvM_err_down2), sf_jet30_csvM[1], sf_jet30_csvM[2], sf_jet30_csvM[3], sf_jet30_csvM[4],
	 sf_jet30_csvM[7], sf_jet30_csvM[8], sf_jet30_csvM[9], sf_jet30_csvM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvT[0], sqrt(sf_jet30_csvT_err_up2), sqrt(sf_jet30_csvT_err_down2), sf_jet30_csvT[1], sf_jet30_csvT[2], sf_jet30_csvT[3], sf_jet30_csvT[4],
	 sf_jet30_csvT[7], sf_jet30_csvT[8], sf_jet30_csvT[9], sf_jet30_csvT[10] );
  }
  else{
  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvL[0], sqrt(sf_jet30_csvL_err_up2), sqrt(sf_jet30_csvL_err_down2), sf_jet30_csvL[1], sf_jet30_csvL[2], sf_jet30_csvL[5], sf_jet30_csvL[6],
	 sf_jet30_csvL[11], sf_jet30_csvL[12], sf_jet30_csvL[13], sf_jet30_csvL[14] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvM[0], sqrt(sf_jet30_csvM_err_up2), sqrt(sf_jet30_csvM_err_down2), sf_jet30_csvM[1], sf_jet30_csvM[2], sf_jet30_csvM[5], sf_jet30_csvM[6],
	 sf_jet30_csvM[11], sf_jet30_csvM[12], sf_jet30_csvM[13], sf_jet30_csvM[14] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30_csvT[0], sqrt(sf_jet30_csvT_err_up2), sqrt(sf_jet30_csvT_err_down2), sf_jet30_csvT[1], sf_jet30_csvT[2], sf_jet30_csvT[5], sf_jet30_csvT[6],
	 sf_jet30_csvT[11], sf_jet30_csvT[12], sf_jet30_csvT[13], sf_jet30_csvT[14] );
  }
  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT > 20 " << std::endl;
  std::cout << " " << std::endl;
  if(isHF){
  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvL[0], sqrt(sf_jet20_csvL_err_up2), sqrt(sf_jet20_csvL_err_down2), sf_jet20_csvL[1], sf_jet20_csvL[2], sf_jet20_csvL[3], sf_jet20_csvL[4],
	 sf_jet20_csvL[7], sf_jet20_csvL[8], sf_jet20_csvL[9], sf_jet20_csvL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvM[0], sqrt(sf_jet20_csvM_err_up2), sqrt(sf_jet20_csvM_err_down2), sf_jet20_csvM[1], sf_jet20_csvM[2], sf_jet20_csvM[3], sf_jet20_csvM[4],
	 sf_jet20_csvM[7], sf_jet20_csvM[8], sf_jet20_csvM[9], sf_jet20_csvM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvT[0], sqrt(sf_jet20_csvT_err_up2), sqrt(sf_jet20_csvT_err_down2), sf_jet20_csvT[1], sf_jet20_csvT[2], sf_jet20_csvT[3], sf_jet20_csvT[4],
	 sf_jet20_csvT[7], sf_jet20_csvT[8], sf_jet20_csvT[9], sf_jet20_csvT[10] );
  }
  else{
  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvL[0], sqrt(sf_jet20_csvL_err_up2), sqrt(sf_jet20_csvL_err_down2), sf_jet20_csvL[1], sf_jet20_csvL[2], sf_jet20_csvL[5], sf_jet20_csvL[6],
	 sf_jet20_csvL[11], sf_jet20_csvL[12], sf_jet20_csvL[13], sf_jet20_csvL[14] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvM[0], sqrt(sf_jet20_csvM_err_up2), sqrt(sf_jet20_csvM_err_down2), sf_jet20_csvM[1], sf_jet20_csvM[2], sf_jet20_csvM[5], sf_jet20_csvM[6],
	 sf_jet20_csvM[11], sf_jet20_csvM[12], sf_jet20_csvM[13], sf_jet20_csvM[14] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : HF (%.4f, %.4f) : LFStats1 (%.4f, %.4f) : LFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet20_csvT[0], sqrt(sf_jet20_csvT_err_up2), sqrt(sf_jet20_csvT_err_down2), sf_jet20_csvT[1], sf_jet20_csvT[2], sf_jet20_csvT[5], sf_jet20_csvT[6],
	 sf_jet20_csvT[11], sf_jet20_csvT[12], sf_jet20_csvT[13], sf_jet20_csvT[14] );
  }
  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  //  if( !isCSV && isHF ){
  if( isHF ){
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT [30, 50) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30to50_cmvaL[0], sqrt(sf_jet30to50_cmvaL_err_up2), sqrt(sf_jet30to50_cmvaL_err_down2), sf_jet30to50_cmvaL[1], sf_jet30to50_cmvaL[2], sf_jet30to50_cmvaL[3], sf_jet30to50_cmvaL[4],
	 sf_jet30to50_cmvaL[7], sf_jet30to50_cmvaL[8], sf_jet30to50_cmvaL[9], sf_jet30to50_cmvaL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30to50_cmvaM[0], sqrt(sf_jet30to50_cmvaM_err_up2), sqrt(sf_jet30to50_cmvaM_err_down2), sf_jet30to50_cmvaM[1], sf_jet30to50_cmvaM[2], sf_jet30to50_cmvaM[3], sf_jet30to50_cmvaM[4],
	 sf_jet30to50_cmvaM[7], sf_jet30to50_cmvaM[8], sf_jet30to50_cmvaM[9], sf_jet30to50_cmvaM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet30to50_cmvaT[0], sqrt(sf_jet30to50_cmvaT_err_up2), sqrt(sf_jet30to50_cmvaT_err_down2), sf_jet30to50_cmvaT[1], sf_jet30to50_cmvaT[2], sf_jet30to50_cmvaT[3], sf_jet30to50_cmvaT[4],
	 sf_jet30to50_cmvaT[7], sf_jet30to50_cmvaT[8], sf_jet30to50_cmvaT[9], sf_jet30to50_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT [50, 70) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet50to70_cmvaL[0], sqrt(sf_jet50to70_cmvaL_err_up2), sqrt(sf_jet50to70_cmvaL_err_down2), sf_jet50to70_cmvaL[1], sf_jet50to70_cmvaL[2], sf_jet50to70_cmvaL[3], sf_jet50to70_cmvaL[4],
	 sf_jet50to70_cmvaL[7], sf_jet50to70_cmvaL[8], sf_jet50to70_cmvaL[9], sf_jet50to70_cmvaL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet50to70_cmvaM[0], sqrt(sf_jet50to70_cmvaM_err_up2), sqrt(sf_jet50to70_cmvaM_err_down2), sf_jet50to70_cmvaM[1], sf_jet50to70_cmvaM[2], sf_jet50to70_cmvaM[3], sf_jet50to70_cmvaM[4],
	 sf_jet50to70_cmvaM[7], sf_jet50to70_cmvaM[8], sf_jet50to70_cmvaM[9], sf_jet50to70_cmvaM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet50to70_cmvaT[0], sqrt(sf_jet50to70_cmvaT_err_up2), sqrt(sf_jet50to70_cmvaT_err_down2), sf_jet50to70_cmvaT[1], sf_jet50to70_cmvaT[2], sf_jet50to70_cmvaT[3], sf_jet50to70_cmvaT[4],
	 sf_jet50to70_cmvaT[7], sf_jet50to70_cmvaT[8], sf_jet50to70_cmvaT[9], sf_jet50to70_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT [70, 100) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet70to100_cmvaL[0], sqrt(sf_jet70to100_cmvaL_err_up2), sqrt(sf_jet70to100_cmvaL_err_down2), sf_jet70to100_cmvaL[1], sf_jet70to100_cmvaL[2], sf_jet70to100_cmvaL[3], sf_jet70to100_cmvaL[4],
	 sf_jet70to100_cmvaL[7], sf_jet70to100_cmvaL[8], sf_jet70to100_cmvaL[9], sf_jet70to100_cmvaL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet70to100_cmvaM[0], sqrt(sf_jet70to100_cmvaM_err_up2), sqrt(sf_jet70to100_cmvaM_err_down2), sf_jet70to100_cmvaM[1], sf_jet70to100_cmvaM[2], sf_jet70to100_cmvaM[3], sf_jet70to100_cmvaM[4],
	 sf_jet70to100_cmvaM[7], sf_jet70to100_cmvaM[8], sf_jet70to100_cmvaM[9], sf_jet70to100_cmvaM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet70to100_cmvaT[0], sqrt(sf_jet70to100_cmvaT_err_up2), sqrt(sf_jet70to100_cmvaT_err_down2), sf_jet70to100_cmvaT[1], sf_jet70to100_cmvaT[2], sf_jet70to100_cmvaT[3], sf_jet70to100_cmvaT[4],
	 sf_jet70to100_cmvaT[7], sf_jet70to100_cmvaT[8], sf_jet70to100_cmvaT[9], sf_jet70to100_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btaggerStr << " "<< flavorStr <<" scale factors for jet pT [100, inf) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sL SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet100toInf_cmvaL[0], sqrt(sf_jet100toInf_cmvaL_err_up2), sqrt(sf_jet100toInf_cmvaL_err_down2), sf_jet100toInf_cmvaL[1], sf_jet100toInf_cmvaL[2], sf_jet100toInf_cmvaL[3], sf_jet100toInf_cmvaL[4],
	 sf_jet100toInf_cmvaL[7], sf_jet100toInf_cmvaL[8], sf_jet100toInf_cmvaL[9], sf_jet100toInf_cmvaL[10] );

  printf("  %sM SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet100toInf_cmvaM[0], sqrt(sf_jet100toInf_cmvaM_err_up2), sqrt(sf_jet100toInf_cmvaM_err_down2), sf_jet100toInf_cmvaM[1], sf_jet100toInf_cmvaM[2], sf_jet100toInf_cmvaM[3], sf_jet100toInf_cmvaM[4],
	 sf_jet100toInf_cmvaM[7], sf_jet100toInf_cmvaM[8], sf_jet100toInf_cmvaM[9], sf_jet100toInf_cmvaM[10] );

  printf("  %sT SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btaggerStr.Data(),
	 sf_jet100toInf_cmvaT[0], sqrt(sf_jet100toInf_cmvaT_err_up2), sqrt(sf_jet100toInf_cmvaT_err_down2), sf_jet100toInf_cmvaT[1], sf_jet100toInf_cmvaT[2], sf_jet100toInf_cmvaT[3], sf_jet100toInf_cmvaT[4],
	 sf_jet100toInf_cmvaT[7], sf_jet100toInf_cmvaT[8], sf_jet100toInf_cmvaT[9], sf_jet100toInf_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  }

      
  std::cout << "Done." << std::endl;

}
