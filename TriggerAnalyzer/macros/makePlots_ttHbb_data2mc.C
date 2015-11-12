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


void makePlots_ttHbb_data2mc( bool printPDF_ = false ){

  TH1::SetDefaultSumw2();

  int NumSamples = 11;
  TFile* file[NumSamples];
  file[0] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_SingleLepton_Run2015D_254231_258158_histo.root");
  file[8] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_EWK_histo.root");
  file[7] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttV_histo.root");
  file[6] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_singlet_histo.root");
  file[5] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttlf_histo.root");
  file[4] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttcc_histo.root");
  file[3] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttb_histo.root");
  file[2] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_tt2b_histo.root");
  file[1] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttbb_histo.root");
  file[9] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttHnonbb_histo.root");
  file[10] = new TFile("HistoFiles/ttHbb_data2mc_treeReader_ttHTobb_histo.root");



  std::vector<TString> histLabels(NumSamples);
  histLabels[0] = "Data";
  histLabels[8] = "EWK";
  histLabels[7] = "tt+W,Z";
  histLabels[6] = "single t";
  histLabels[5] = "tt+lf";
  histLabels[4] = "tt+cc";
  histLabels[3] = "tt+b";
  histLabels[2] = "tt+2b";
  histLabels[1] = "tt+bb";
  histLabels[9] = "ttHnonbb";
  histLabels[10] = "ttHbb";


  Color_t color[NumSamples];
  color[0] = kBlack;
  color[8] = kAzure+2;
  color[7] = kBlue-10;
  color[6] = kMagenta;
  color[5] = kRed-7;
  color[4] = kRed+1;
  color[3] = kRed-2;
  color[2] = kRed+0;
  color[1] = kRed+3;
  color[9] = kGreen+2;
  color[10] = kBlue;




  TString dirprefix = "Images/Images_2015_10_20_ttHbb_data2mc/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

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



  std::vector<std::string> histoname1;

  histoname1.push_back("h_lepton_pt");
  histoname1.push_back("h_lepton_eta");
  histoname1.push_back("h_lepton_phi");

  histoname1.push_back("h_jet_pt");
  histoname1.push_back("h_jet_eta");
  histoname1.push_back("h_jet_phi");
  histoname1.push_back("h_jet_csv");

  histoname1.push_back("h_met_pt");
  histoname1.push_back("h_met_phi");
  histoname1.push_back("h_mht_pt");
  histoname1.push_back("h_mht_phi");

  histoname1.push_back("h_HT");


  std::vector<std::string> histoname2;
  histoname2.push_back("h_category_yield");
  histoname2.push_back("h_category_yield_1e");
  histoname2.push_back("h_category_yield_1m");

  histoname2.push_back("h_event_selection");
  histoname2.push_back("h_mu_event_selection");
  histoname2.push_back("h_ele_event_selection");

  histoname2.push_back("h_numJet");
  histoname2.push_back("h_numBtag");

  histoname2.push_back("h_deltaR_jet_lep");


  for( int iHist=0; iHist<int(histoname1.size()); iHist++ ){
    for( int iCat=0; iCat<NumCat; iCat++ ){
      std::string suffix = "_" + cat_labels[iCat];
      std::string new_name = histoname1[iHist] + suffix;
      histoname2.push_back(new_name);
    }
  }


  std::vector<std::string> histoname3;
  histoname3.push_back("h_category_yield");
  histoname3.push_back("h_category_yield_1e");
  histoname3.push_back("h_category_yield_1m");


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGaxis::SetMaxDigits(4);

  TString lumiinfo = "553 pb^{-1} (13 TeV)";
  TLatex LumiInfoLatex(0.70, 0.91, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.04);

  //TString cmsinfo =   "CMS Preliminary";
  TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex(0.185, 0.91, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextFont(61);
  CMSInfoLatex.SetTextSize(0.055); //SBOUTLE

  std::string publishinfo =   "Preliminary"; //DPUIGH
  TLatex PublishInfoLatex(0.285, 0.91, publishinfo.c_str()); //SBOUTLE
  PublishInfoLatex.SetNDC();
  PublishInfoLatex.SetTextFont(52);
  PublishInfoLatex.SetTextSize(0.045); //SBOUTLE



  double scale_ttH = 30.;


  TString plotname;


  TCanvas* myC1 = new TCanvas("myC1", "myC1", 600,700);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1.e-5;
  myC1->Divide(1,2,small,small);
  const float padding=1e-5; const float ydivide=0.3;
  myC1->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
  myC1->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
  // myC1->GetPad(1)->SetLeftMargin(.11);
  // myC1->GetPad(2)->SetLeftMargin(.11);
  myC1->GetPad(1)->SetLeftMargin(.12);
  myC1->GetPad(2)->SetLeftMargin(.12);
  myC1->GetPad(1)->SetRightMargin(.05);
  myC1->GetPad(2)->SetRightMargin(.05);
  myC1->GetPad(1)->SetBottomMargin(.3);
  myC1->GetPad(2)->SetBottomMargin(.3);
  myC1->GetPad(1)->Modified();
  myC1->GetPad(2)->Modified();
  myC1->cd(1);
  gPad->SetBottomMargin(small);
  gPad->Modified();



  for( int i=0; i<int(histoname2.size()); i++ ){

    TString temp = histoname2[i];
    
    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    // //TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
    // TLegend *legend = new TLegend(0.13,0.83,0.88,0.89);

    // legend->SetFillColor(kWhite);
    // legend->SetLineColor(kWhite);
    // legend->SetShadowColor(kWhite);
    // legend->SetTextFont(42);
    // legend->SetTextSize(0.04);

    // legend->SetNColumns(4);

    TLegend *legend = new TLegend(0.70,0.50,0.84,0.89);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);


    TLegend *legend_left = new TLegend(0.15,0.50,0.32,0.89);

    legend_left->SetFillColor(kWhite);
    legend_left->SetLineColor(kWhite);
    legend_left->SetShadowColor(kWhite);
    legend_left->SetTextFont(42);
    legend_left->SetTextSize(0.04);



    int rebin = 1;

    // if( temp.Contains("_mass") ) rebin = 2;
    // if( temp.Contains("_L1HTT") ) rebin = 2;
    // if( temp.Contains("_met") ) rebin = 5;
    // if( temp.Contains("_jet_") && temp.Contains("_pt") ) rebin = 10;
    // if( temp.Contains("_jet_") && temp.Contains("_eta") ) rebin = 4;
    // if( temp.Contains("_jet_") && temp.Contains("_phi") ) rebin = 4;
    // if( temp.Contains("_csv") ) rebin = 8;

    // if( temp.Contains("_diele_mass_closestZmass") ) rebin = 5;

    if( temp.Contains("_HT") && !temp.Contains("4j") ) rebin = 2;

    TH1D* hist[NumSamples];
    TH1D* hist_sum = NULL;
    bool firstFill = true;
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      hist[iSample] = (TH1D*)file[iSample]->Get(temp);
      hist[iSample]->Rebin(rebin);

      if( iSample>0 ){
	hist[iSample]->SetLineColor(color[iSample]);
	if( iSample<NumSamples-2 ) hist[iSample]->SetFillColor(color[iSample]);
	else hist[iSample]->SetLineWidth(2);

	if( firstFill && iSample>0 ){
	  firstFill = false;
	  hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp);
	}
	else if( iSample<NumSamples-2 ){
	  hist_sum->Add( hist[iSample] );
	}
      }
      else {
	hist[iSample]->SetMarkerStyle(20);
	//hist[iSample]->SetLineWidth(20);
      }

      if( iSample==NumSamples-1 || iSample==NumSamples-2 ) hist[iSample]->Scale(scale_ttH);

      if( iSample==0 )                legend->AddEntry(hist[iSample],histLabels[iSample],"pe1");
      else if( iSample<NumSamples-2 ) legend->AddEntry(hist[iSample],histLabels[iSample],"f");
      else                            legend->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");

      if( iSample==0 )                legend_left->AddEntry(hist[iSample],histLabels[iSample],"pe1");
      else if( iSample<NumSamples-2 ) legend_left->AddEntry(hist[iSample],histLabels[iSample],"f");
      else                            legend_left->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");
    }// end loop over samples


    THStack *hs = new THStack("hs","");
    for( int iSample=NumSamples-1; iSample>-1; iSample-- ){

      if( iSample>0 && iSample<NumSamples-2 ) hs->Add(hist[iSample]);
    }


    // double ratioMax = 1.6;
    // double ratioMin = 0.5;
    double ratioMax = 2.3;
    double ratioMin = 0.0;

    int nbins = hist[0]->GetNbinsX();

    double xmin = hist[0]->GetBinLowEdge(1);
    double xmax = hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins);

    TH1D* myRatio = new TH1D("ratio", "", nbins, xmin, xmax );

    myRatio->SetStats(0);
    myRatio->Sumw2();
    myRatio->SetLineColor(kBlack);
    myRatio->SetMarkerColor(kBlack);
    myRatio->Divide(hist[0],hist_sum);

    myRatio->SetMinimum(ratioMin);
    myRatio->SetMaximum(ratioMax);
    // double ratioMax = 2.3;
    // double ratioMin = 0.0;
    myRatio->GetYaxis()->SetNdivisions(50000+404);
    // double ratioMax = 1.6;
    // double ratioMin = 0.5;
    //myRatio->GetYaxis()->SetNdivisions(50000+204);
    myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetTitleOffset(1.1);
    myRatio->GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle()); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.12);
    myRatio->GetXaxis()->SetLabelOffset(0.04);
    myRatio->GetXaxis()->SetTitleSize(0.12);
    myRatio->GetYaxis()->SetTitle("Data/Bkg");
    myRatio->GetYaxis()->SetTitleSize(0.09);
    myRatio->GetYaxis()->SetTitleOffset(.55);
    myC1->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();

    myRatio->GetYaxis()->CenterTitle(kTRUE);


    if( temp.Contains("_category_yield") || temp.Contains("_selection") ){
      for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,hist[0]->GetXaxis()->GetBinLabel(iBin+1));
    }


    hist[0]->SetStats(0);

    hist[0]->GetYaxis()->SetTitleOffset(1.2);
    hist[0]->GetYaxis()->SetTitleSize(0.05);

    hist[0]->GetYaxis()->SetTitle("Number of Events");


    int max_bin_data = hist[0]->GetMaximumBin();
    double max_data = hist[0]->GetBinContent(max_bin_data) + hist[0]->GetBinError(max_bin_data);

    int max_bin_mc = hist_sum->GetMaximumBin();
    double max_mc = hist_sum->GetBinContent(max_bin_mc) + hist_sum->GetBinError(max_bin_mc);

    double max_content = std::max(max_data, max_mc);

    hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * max_content);


    bool rerange = false;

    if( (temp.Contains("_jet_csv")) ){ rerange=true; xmin = -0.04; xmax = 1.048; }
    if( (temp.Contains("_HT")) ){ rerange=true; xmin = 100+0.0001; xmax = 1400-0.001; }
    if( (temp.Contains("_mht_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 200-0.001; }
    if( (temp.Contains("_met_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 400-0.001; }
    //if( (temp.Contains("_category")) ){ rerange=true; xmin = 0.001; xmax = 6.99; }

    if( rerange ){
      hist[0]->GetXaxis()->SetRangeUser(xmin,xmax);
      myRatio->GetXaxis()->SetRangeUser(xmin,xmax);
    }

    // if( temp.Contains("_mass") ){
    //   hist[0]->GetXaxis()->SetRangeUser(40.,140.);
    //   myRatio->GetXaxis()->SetRangeUser(40.,140.);
    // }

    // if( temp.Contains("_met") ){
    //   hist[0]->GetXaxis()->SetRangeUser(0.,150.);
    //   myRatio->GetXaxis()->SetRangeUser(0.,150.);
    // }
    // // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);


    TLine* myLine;
    myLine = new TLine(xmin, 1, xmax, 1);


    // Plot
    myC1->cd(1);
    hist[0]->Draw("pe1");
    hs->Draw("histsame");
    hist[NumSamples-2]->Draw("histsame");
    hist[NumSamples-1]->Draw("histsame");
    hist[0]->Draw("pe1same");

    if( temp.Contains("_jet_csv") ) legend_left->Draw();
    else                            legend->Draw();

    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    myC1->cd(2);
    myRatio->SetLineWidth(2);
    myRatio->Draw("pe1");
    myLine->Draw();

    myC1->GetPad(1)->SetLogy(0);

    myC1->GetPad(1)->RedrawAxis();
    myC1->GetPad(2)->RedrawAxis();

    plotname = dirprefix + temp_mh + "_data2mc_lin.png";
    myC1->Print(plotname);

    plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
    if( printPDF_ ) myC1->Print(plotname);


    if( temp.Contains("_category_yield") || temp.Contains("h_num") ){
      // log
      hist[0]->GetYaxis()->SetRangeUser(0.4,12 * max_content);
      if( temp.Contains("_selection") ){
	hist[0]->GetYaxis()->SetRangeUser(1000,12 * max_content);
      }

      myC1->cd(1);
      //gPad->SetLogy(1);

      myC1->GetPad(1)->SetLogy(1);

      myC1->cd(1);
      hist[0]->Draw("pe1");
      hs->Draw("histsame");
      hist[NumSamples-2]->Draw("histsame");
      hist[NumSamples-1]->Draw("histsame");
      hist[0]->Draw("pe1same");

      if( temp.Contains("_jet_csv") ) legend_left->Draw();
      else                            legend->Draw();

      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("pe1");
      myLine->Draw();

      myC1->GetPad(1)->RedrawAxis();
      myC1->GetPad(2)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_log.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
      if( printPDF_ ) myC1->Print(plotname);

      //gPad->SetLogy(0);
    }


    delete myRatio;
    delete myLine;
    delete legend;
  } // end loop on hists




  for( int i=0; i<int(histoname3.size()); i++ ){

    TString temp = histoname3[i];
    
    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    TH1D* hist[NumSamples];
    TH1D* hist_sum = NULL;
    bool firstFill = true;
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      hist[iSample] = (TH1D*)file[iSample]->Get(temp);

      if( iSample>0 ){
	if( firstFill && iSample>0 ){
	  firstFill = false;
	  hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp);
	}
	else if( iSample<NumSamples-2 ){
	  hist_sum->Add( hist[iSample] );
	}
      }
    }// end loop over samples

    int NumBins = hist[0]->GetNbinsX();

    double xmin = hist[0]->GetBinLowEdge(1);
    double xmax = hist[0]->GetBinLowEdge(NumBins) + hist[0]->GetBinWidth(NumBins);

    TH1D* myRatio = new TH1D("ratio_"+temp, "", NumBins, xmin, xmax );
    myRatio->Divide(hist[0],hist_sum);


    std::cout << " ============================================= " << std::endl;

    std::cout << " \t " << temp_mh << std::endl;

    std::cout << "\\begin{tabular}{|l|c|c|c|c|c|c|c|c|c|c|} \\hline" << std::endl; 
    std::cout << "\t";
    for( int iBin=0; iBin<NumBins; iBin++ ){
      std::cout << "&\t" << hist[0]->GetXaxis()->GetBinLabel(iBin+1);
    }
    std::cout << "\\\\ \\hline \\hline" << std::endl;

    printf("%8s", histLabels[NumSamples-1].Data());
    //std::cout << histLabels[NumSamples-1];
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist[NumSamples-1]->GetBinContent(iBin+1)/scale_ttH, hist[NumSamples-1]->GetBinError(iBin+1)/scale_ttH);
    }
    std::cout << "\\\\" << std::endl;    

    printf("%8s", histLabels[NumSamples-2].Data());
    //std::cout << histLabels[NumSamples-2];
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist[NumSamples-2]->GetBinContent(iBin+1)/scale_ttH, hist[NumSamples-2]->GetBinError(iBin+1)/scale_ttH);
    }
    std::cout << " " << std::endl;  
    std::cout << "\\\\ \\hline" << std::endl;

    for( int iSample=1; iSample<NumSamples-2; iSample++ ){
      //std::cout << histLabels[iSample];
      printf("%8s", histLabels[iSample].Data());
      for( int iBin=0; iBin<NumBins; iBin++ ){
	printf(" & %6.1f $\\pm$ %3.1f", hist[iSample]->GetBinContent(iBin+1), hist[iSample]->GetBinError(iBin+1));
      }
      std::cout << " \\\\" << std::endl; 
    }
    std::cout << "\\hline" << std::endl;

    //std::cout << "Data";
    printf("%8s", "Data");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.0f ", hist[0]->GetBinContent(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  

    printf("%8s", "TotBkg");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist_sum->GetBinContent(iBin+1), hist_sum->GetBinError(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  


    //std::cout << "Data/Bkg";
    printf("%8s", "Data/Bkg");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.2f $\\pm$ %3.2f", myRatio->GetBinContent(iBin+1), myRatio->GetBinError(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  

    std::cout << " \\end{tabular} " << std::endl;


    delete myRatio;
  }
  std::cout << " ============================================= " << std::endl;





  // Close the files
  for( int iFile=0; iFile<NumSamples; iFile++ ) file[iFile]->Close();
  std::cout << " Done! " << std::endl;
}
