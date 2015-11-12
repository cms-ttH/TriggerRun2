#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TKey.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>

#include "TMath.h"

//*****************************************************************************

//*****************************************************************************


void calculate_PU_reweight( TString dataFileName = "", TString mcFileName = "" ){

  TH1::SetDefaultSumw2();

  TFile* file_data = new TFile(dataFileName);
  TFile* file_mc   = new TFile(mcFileName);

  TH1D* h_data = (TH1D*)file_data->Get("h_numPVs");
  TH1D* h_mc = (TH1D*)file_mc->Get("h_numPVs");
 
  h_data->Scale(1. / h_data->Integral());
  h_mc->Scale(1. / h_mc->Integral());
 
  TH1D* h_ratio = (TH1D*)h_data->Clone("h_ratio");
  h_ratio->Divide(h_mc);
 
  for( int iBin=0; iBin<h_ratio->GetNbinsX(); iBin++ ) std::cout << "  PUscale[" << iBin << "] = " << h_ratio->GetBinContent(iBin+1) << ";" << std::endl;
 

  // Close the file
  file_data->Close();
  file_mc->Close();
}
