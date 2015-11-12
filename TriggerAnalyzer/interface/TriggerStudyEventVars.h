#ifndef TriggerRun2_TriggerAnalyzer_TriggerStudyEventVars_h
#define TriggerRun2_TriggerAnalyzer_TriggerStudyEventVars_h

//
// Dependencies (#includes)
//
#include <iostream>
#include <vector>
#include "TLorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

using namespace std;



typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<double> vdouble;
typedef std::vector<string> vstring;
typedef std::vector<bool> vbool;
typedef std::vector<int> vint;
typedef std::vector< TLorentzVector > vecTLorentzVector;

//
// Utility Class for Handling Event Variables
//

struct triggerStudyEventVars{


  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  explicit triggerStudyEventVars() { }

  int run_;
  int lumi_;
  long evt_;

  double  wgt_generator_;

  int numJets_;
  int numTags_;
  int numPVs_;
  int numTruePVs_;
  int numGenPVs_;

  int additionalJetEventId_;

  /////

  int pass_HLT_Ele27_eta2p1_WP75_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v_;
  int pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v_;

  int pass_HLT_Ele25WP60_SC4_Mass55_v_;
  int pass_HLT_Ele27_WP85_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_;

  int pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v_;
  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v_;

  int pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;

  /////


  double pfMET_pt_;
  double pfMET_phi_;
  double L1HTT_;
  double L1HTT_bxm2_;
  double L1HTT_bxm1_;
  double L1HTT_bxp1_;
  double L1HTT_bxp2_;

  vdouble jet_pt_;
  vdouble jet_eta_;
  vdouble jet_phi_;
  vdouble jet_energy_;
  vdouble jet_csv_;
  vint    jet_flavor_;

  vdouble jet_nocc_pt_;
  vdouble jet_nocc_eta_;
  vdouble jet_nocc_phi_;
  vdouble jet_nocc_energy_;
  vdouble jet_nocc_csv_;
  vint    jet_nocc_flavor_;

  double mass_leplep_;
  double dR_leplep_;
  Int_t  oppositeLepCharge_;

  vint lepton_trkCharge_;
  vint lepton_isMuon_;
  vint lepton_isTight_;
  vint lepton_isLoose_;
  vint lepton_isPhys14L_;
  vint lepton_isPhys14M_;
  vint lepton_isPhys14T_;
  vint lepton_isSpring15L_;
  vint lepton_isSpring15M_;
  vint lepton_isSpring15T_;
  vint lepton_genId_;
  vint lepton_genParentId_;
  vint lepton_genGrandParentId_;
  vdouble lepton_pt_;
  vdouble lepton_eta_;
  vdouble lepton_phi_;
  vdouble lepton_energy_;
  vdouble lepton_relIso_;
  vdouble lepton_relIsoR04_;
  vdouble lepton_iso_sumChargedHadronPt_;
  vdouble lepton_iso_sumNeutralHadronEt_;
  vdouble lepton_iso_sumPhotonEt_;
  vdouble lepton_iso_sumPUPt_;
  vdouble lepton_mvaTrigValue_;
  vdouble lepton_scSigmaIEtaIEta_;
  vdouble lepton_full5x5_scSigmaIEtaIEta_;
  vdouble lepton_hadronicOverEm_;
  vdouble lepton_relEcalIso_;
  vdouble lepton_relHcalIso_;
  vdouble lepton_relTrackIso_;
  vdouble lepton_OneOESuperMinusOneOP_;
  vint lepton_numMissingHits_;
  vint lepton_isEB_;
  vint lepton_passHLTId_;
  vint lepton_passConversionVeto_;
  vint lepton_inCrack_;
  vdouble lepton_scEta_;
  vdouble lepton_dEtaSCTrackAtVtx_;
  vdouble lepton_dPhiSCTrackAtVtx_;
  vdouble lepton_d0_;
  vdouble lepton_dZ_;
  vint lepton_isGlobalMuon_;
  vint lepton_isTrackerMuon_;
  vint lepton_isPFMuon_;
  vdouble lepton_normalizedChi2_;
  vint lepton_numberOfValidMuonHits_;
  vint lepton_numberOfValidPixelHits_;
  vint lepton_trackerLayersWithMeasurement_;
  vint lepton_numberOfMatchedStations_;

  vint lepton_ele_matchHLT_hltL1sL1SingleEG25_;
  vint lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_;

  vint lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25_;
  vint lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_;

  vint lepton_ele_matchHLT_hltL1sL1EG25erHTT125_;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter_;
  vint lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter_;

  vdouble hltL1SingleEG25_pt_;
  vdouble hltL1SingleEG25_eta_;
  vdouble hltL1SingleEG25_phi_;

  vdouble hltEle27WP85Gsf_pt_;
  vdouble hltEle27WP85Gsf_eta_;
  vdouble hltEle27WP85Gsf_phi_;
  vstring hltEle27WP85Gsf_filter_;

  vdouble hltJet30_pt_;
  vdouble hltJet30_eta_;
  vdouble hltJet30_phi_;

  vdouble hltBtagJet30_pt_;
  vdouble hltBtagJet30_eta_;
  vdouble hltBtagJet30_phi_;

  vdouble hltPFHT200Jet30_pt_;
  vdouble hltPFHT200Jet30_eta_;
  vdouble hltPFHT200Jet30_phi_;
  vint    hltPFHT200Jet30_id_;

  vint    hlt_accept_;
  vstring hlt_name_;

  vint    l1t_accept_;
  vstring l1t_name_;

  void initialize();

};


typedef std::vector<triggerStudyEventVars> vtriggerStudyEventVars;


void triggerStudyEventVars::initialize(){

  run_  = -99;
  lumi_ = -99;
  evt_ = -99;

  wgt_generator_ = -99.9;

  numJets_ = -99;
  numTags_ = -99;
  numPVs_  = -99;
  numTruePVs_ = -99;
  numGenPVs_ = -99;

  additionalJetEventId_ = -99;

  pass_HLT_Ele27_eta2p1_WP75_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v_ = -99;

  pass_HLT_Ele25WP60_SC4_Mass55_v_ = -99;
  pass_HLT_Ele27_WP85_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_ = -99;


  pass_HLT_Ele27_eta2p1_WPTight_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v_ = -99;
  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v_ = -99;

  pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_ = -99;


  pfMET_pt_ = -1;
  pfMET_phi_ = -1;

  L1HTT_ = -1;
  L1HTT_bxm2_ = -1;
  L1HTT_bxm1_ = -1;
  L1HTT_bxp1_ = -1;
  L1HTT_bxp2_ = -1;


  jet_pt_.clear();
  jet_eta_.clear();
  jet_phi_.clear();
  jet_energy_.clear();
  jet_csv_.clear();
  jet_flavor_.clear();

  jet_nocc_pt_.clear();
  jet_nocc_eta_.clear();
  jet_nocc_phi_.clear();
  jet_nocc_energy_.clear();
  jet_nocc_csv_.clear();
  jet_nocc_flavor_.clear();

  mass_leplep_ = -99;
  dR_leplep_ = -99;
  oppositeLepCharge_ = -99;


  lepton_trkCharge_.clear();
  lepton_isMuon_.clear();
  lepton_isTight_.clear();
  lepton_isLoose_.clear();
  lepton_isPhys14L_.clear();
  lepton_isPhys14M_.clear();
  lepton_isPhys14T_.clear();
  lepton_isSpring15L_.clear();
  lepton_isSpring15M_.clear();
  lepton_isSpring15T_.clear();
  lepton_genId_.clear();
  lepton_genParentId_.clear();
  lepton_genGrandParentId_.clear();
  lepton_pt_.clear();
  lepton_eta_.clear();
  lepton_phi_.clear();
  lepton_energy_.clear();
  lepton_relIso_.clear();
  lepton_relIsoR04_.clear();
  lepton_iso_sumChargedHadronPt_.clear();
  lepton_iso_sumNeutralHadronEt_.clear();
  lepton_iso_sumPhotonEt_.clear();
  lepton_iso_sumPUPt_.clear();
  lepton_mvaTrigValue_.clear();
  lepton_scSigmaIEtaIEta_.clear();
  lepton_full5x5_scSigmaIEtaIEta_.clear();
  lepton_hadronicOverEm_.clear();
  lepton_relEcalIso_.clear();
  lepton_relHcalIso_.clear();
  lepton_relTrackIso_.clear();
  lepton_OneOESuperMinusOneOP_.clear();
  lepton_numMissingHits_.clear();
  lepton_isEB_.clear();
  lepton_passHLTId_.clear();
  lepton_passConversionVeto_.clear();
  lepton_inCrack_.clear();
  lepton_scEta_.clear();
  lepton_dEtaSCTrackAtVtx_.clear();
  lepton_dPhiSCTrackAtVtx_.clear();
  lepton_d0_.clear();
  lepton_dZ_.clear();
  lepton_isGlobalMuon_.clear();
  lepton_isTrackerMuon_.clear();
  lepton_isPFMuon_.clear();
  lepton_normalizedChi2_.clear();
  lepton_numberOfValidMuonHits_.clear();
  lepton_numberOfValidPixelHits_.clear();
  lepton_trackerLayersWithMeasurement_.clear();
  lepton_numberOfMatchedStations_.clear();


  lepton_ele_matchHLT_hltL1sL1SingleEG25_.clear();
  lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_.clear();

  lepton_ele_matchHLT_hltL1sL1SingleIsoEG22erOrSingleEG25_.clear();
  lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_.clear();

  lepton_ele_matchHLT_hltL1sL1EG25erHTT125_.clear();
  lepton_ele_matchHLT_hltL1EGHttEle27WP85GsfTrackIsoFilter_.clear();
  lepton_ele_matchHLT_hltL1EGHttEle27WPLooseGsfTrackIsoFilter_.clear();


  hltL1SingleEG25_pt_.clear();
  hltL1SingleEG25_eta_.clear();
  hltL1SingleEG25_phi_.clear();

  hltEle27WP85Gsf_pt_.clear();
  hltEle27WP85Gsf_eta_.clear();
  hltEle27WP85Gsf_phi_.clear();
  hltEle27WP85Gsf_filter_.clear();

  hltJet30_pt_.clear();
  hltJet30_eta_.clear();
  hltJet30_phi_.clear();

  hltBtagJet30_pt_.clear();
  hltBtagJet30_eta_.clear();
  hltBtagJet30_phi_.clear();

  hltPFHT200Jet30_pt_.clear();
  hltPFHT200Jet30_eta_.clear();
  hltPFHT200Jet30_phi_.clear();
  hltPFHT200Jet30_id_.clear();

  hlt_accept_.clear();
  hlt_name_.clear();

  l1t_accept_.clear();
  l1t_name_.clear();


  return;
}

  

#endif
