#ifndef TAUMONTOOL_H
#define TAUMONTOOL_H

//*********************************************************************
// tauMonTool.h : Implementation of offline tau Data Quality 
//                histogramming, inheriting from ManagedMonitorToolBase
//
// authors:  C Cuenca Almenar, S Demers, E Ideal, A Leister, YALE
//*********************************************************************

#include "TH1.h"
#include "TH2.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "tauEvent/TauJetContainer.h"

class tauMonTool : public ManagedMonitorToolBase {
  
 public:
  // constructor
  tauMonTool( const std::string & type, 
	      const std::string & name, 
	      const IInterface* parent );
  
  // destructor
  virtual ~tauMonTool();		
  
  StatusCode bookHistograms();
  StatusCode fillHistograms();
  StatusCode procHistograms();
    StatusCode Book1DHist(TH1 ** hist, MonGroup * monName, std::string histName, std::string histTitle, int NBins, double lowBin, double highBin);
  StatusCode Book2DHist(TH2 ** hist, MonGroup * monName, std::string histName, std::string histTitle, int NXBins, double lowXBin, double highXBin, int NYBins, double lowYBin, double highYBin);
  
 protected:
  
  std::string m_tauJetKey;
  std::string m_tauCommonDetailsContainer;
  std::string m_tauCommonExtraDetailsContainer;
  bool        m_doTrigger;
  
  int  m_maxNLB; 
  int  m_currentLB;
  
  
  //global plots
  struct s_basicPlots {
    TH1* h_ntaus_vs_LB[6];
    int  h_ntausLB[6];
    TH1* h_ntaus;
    
    TH1* h_author; //not used in 2012
    TH1* h_eta;
    TH1* h_phi;
    TH1* h_et;
    TH1* h_charge;
    TH1* h_numTracks;
    TH1* h_nHighPTtaus;
    
    TH2* h_EtVsEta;
    TH2* h_EtVsPhi;
    TH2* h_PhiVsEta;

    TH2* h_Eta_vs_LB;
    TH2* h_Phi_vs_LB;
  } m_basicPlots, m_emPlots, m_jetPlots, m_tauPlots, m_mbtsPlots, m_AtlasReady;
  
  
  //trigger
  TH1* m_triggers;
  
  struct s_BDTFolder{
    TH1* h_et;
    TH1* h_eta;
    TH1* h_phi;
    TH1* h_nTracks;
  };
  
  struct s_idFolder {
    TH1* h_tauCutLoose;
    TH1* h_tauCutMedium;
    TH1* h_tauCutTight;
    TH1* h_tauLlhLoose;
    TH1* h_tauLlhMedium;
    TH1* h_tauLlhTight;
    TH1* h_electronVetoLoose;
    TH1* h_electronVetoMedium;
    TH1* h_electronVetoTight;
    TH1* h_muonVeto;
    TH1* h_eleBDTLoose;
    TH1* h_eleBDTMedium;
    TH1* h_eleBDTTight;    
    TH1* h_BDTJetScore;
    TH1* h_BDTEleScore; 
    TH1* h_BDTJetScoreSigTrans;
    TH1* h_BDTJetScoreBkgTrans;

    s_BDTFolder BDTLooseFolder;
    s_BDTFolder BDTMedFolder;
  };
  
  
  struct s_trkFolder {
    TH1* h_leadTrkPt;
    TH1* h_massTrkSys; 
    TH1* h_trkWidth2; 
    TH1* h_trFlightPathSig;
    TH1* h_ipSigLeadTrk;
    TH1* h_ipZ0SinThetaSigLeadTrk;
    TH1* h_d0;
    TH1* h_z0;
    TH1* h_phi;
    TH1* h_eta;
    TH1* h_pT;
    TH1* h_nHighPTtaus;
    TH1* h_numberOfTRTHighThresholdHits;
    TH1* h_numberOfTRTHits;
    TH1* h_numberOfTRTHighThresholdOutliers;
    TH1* h_numberOfTRTOutliers;
    TH1* h_numberOfSCTHits;
    TH1* h_numberOfPixelHits;	

    TH2* h_z0_vs_LB;
  };
  
  
  struct s_caloFolder {
    TH1* h_eta;
    TH1* h_phi;
    TH1* h_etEMCalib;
    TH1* h_etHadCalib;
    TH1* h_etEMAtEMScale;
    TH1* h_etHadAtEMScale;
    TH1* h_EMRadius;
    TH1* h_centFrac;
    TH1* h_hadRadius;
    TH1* h_isolFrac;
    TH1* h_stripWidth2;
    TH1* h_nStrip;
    TH1* h_trkAvgDist;
    TH1* h_dRmax;
    //TH1* h_cellBasedEnergyRing[7];

    TH2* h_centFrac_vs_LB;
    TH2* h_isolFrac_vs_LB;
  };
  
  
  struct s_combFolder {
    TH1* h_etOverPtLeadTrack;
  };
  
  
  struct s_kinFolder {
    TH1* h_ntaus;  
    TH1* h_author; //not used in 2012
    TH1* h_eta;
    TH1* h_phi;
    TH1* h_et;
    TH1* h_charge;

    TH2* h_PhiVsEta;
    TH2* h_Eta_vs_LB;
    TH2* h_Phi_vs_LB;
  };
  
  struct s_mainFolder {
    s_idFolder   idFolder;
    s_trkFolder  trkFolder;	  
    s_caloFolder caloFolder;
    s_combFolder combFolder;
    s_kinFolder  kinFolder;
    
  } m_tauB, m_tauCR, m_tauE;
  
  
  // Z details
  TH1* m_eta_Tau_Z;
  TH2* m_pTVsEta_Tau_Z;
  TH2* m_pTVsEta_Lepton_Z;

  // W details
  TH2* m_pTVsEta_Tau_W;
  
  // Merged Candidates
    /*
  TH1* m_etaCalo_MERG;
  TH1* m_phiCalo_MERG;
  TH1* m_emRadius_MERG;
  TH1* m_isolationFraction_MERG;
  TH1* m_centralityFraction_MERG;
  TH1* m_stripWidth2_MERG;
  TH1* m_numStripCells_MERG;
  TH1* m_etEM_MERG;
  TH1* m_etEMCalib_MERG;
  TH1* m_etHad_MERG;
  TH1* m_etHadCalib_MERG;
  TH1* m_etOverPtLeadTrack_MERG;
  TH1* m_leadingTrackPT_MERG;
  TH1* m_trFlightPathSig_MERG;
  TH1* m_ipSigLeadTrack_MERG;
  TH1* m_tauCutLoose_MERG;
  TH1* m_tauCutMedium_MERG;
  TH1* m_tauCutTight_MERG;
  TH1* m_electronVeto_MERG;
  TH1* m_muonVeto_MERG;
  TH1* m_tauLlhTight_MERG;
  TH1* m_tauLlhMedium_MERG;
  TH1* m_tauLlhLoose_MERG;
  TH2* m_EtOverPtLTVsEta_MERG;
  TH2* m_EtOverPtLTVsPhi_MERG;
  TH2* m_etEmVsEta_MERG;
  TH2* m_etHadVsEta_MERG;
  TH2* m_etEmVsPhi_MERG;
  TH2* m_etHadVsPhi_MERG;
  TH1* m_etChrgHAD_MERG;
  TH1* m_etIsolEM_MERG;
  TH1* m_etIsolHAD_MERG;
  TH1* m_nAssocTracksCore_MERG;
  TH1* m_nAssocTracksIsol_MERG;
  TH1* m_massTrk3P_MERG;
  TH1* m_rWidth2Trk3P_MERG;
  TH1* m_signD0Trk3P_MERG;
  TH1* m_etHadAtEMScale_MERG;
  TH1* m_etEMAtEMScale_MERG;
  TH1* m_etEMCL_MERG;
  TH1* m_etChrgEM_MERG;
  TH1* m_etNeuEM_MERG;
  TH1* m_etResNeuEM_MERG;
  TH1* m_numPi0_MERG;
  TH1* m_etChrgHADoverPttot_MERG;
  TH1* m_z0SinThetaSig_MERG;
  TH1* m_etIsolFrac_MERG;
  TH1* m_HighTRTHits_MERG;
  TH1* m_LowTRTHits_MERG;
  TH1* m_SCTHits_MERG;
  TH1* m_PixelHits_MERG;
  TH2* m_numPi0VsEta_MERG;
  TH2* m_numPi0VsPhi_MERG;
  TH2* m_nTracksCoreVsEta_MERG;
  TH2* m_nTracksCoreVsPhi_MERG;
  TH2* m_etIsolFracVsEta_MERG;
  TH2* m_etIsolFracVsPhi_MERG;
  */
 private:

  StoreGateSvc* m_storeGate;
  ToolHandle<Trig::TrigDecisionTool> m_trigDec;
	
  StatusCode bookBasicPlots(s_basicPlots& trigPlots, MonGroup &aGroup, std::string preffix);
  void fillBasicPlots(s_basicPlots& someBasicPlots, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau);
  
  StatusCode bookHistos(s_mainFolder& mainFolder, std::string folderName, Interval_t interval);
  StatusCode bookKinHistos(s_kinFolder& folder,  MonGroup &aGroup);
  StatusCode bookIDHistos(s_idFolder& folder,std::string folderName, Interval_t interval) ;
  StatusCode bookTrackHistos(s_trkFolder& folder,std::string folderName, Interval_t interval);
  StatusCode bookCaloHistos(s_caloFolder& folder,std::string folderName, Interval_t interval);
  StatusCode bookCombHistos(s_combFolder& folder,std::string folderName, Interval_t interval);
  StatusCode bookBDTLooseHistos(s_BDTFolder& folder,std::string folderName, Interval_t interval);
  StatusCode bookBDTMedHistos(s_BDTFolder& folder,std::string folderName, Interval_t interval);
  StatusCode bookPhysicsHistograms();
  //StatusCode bookOldHistograms();
  
  StatusCode fillHistograms(s_mainFolder& mainFolder, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau, 
			    const Analysis::TauCommonDetails* p_taudetails,const Analysis::TauPID* p_tauid, const Trk::TrackSummary* trackSummary,
			    const Trk::Perigee* perigee);
  StatusCode fillKinHistos(s_kinFolder& folder, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau) ;
  StatusCode fillIDHistos(s_idFolder& folder, const Analysis::TauPID* p_tauid, Analysis::TauJetContainer::const_iterator ftau);
  StatusCode fillTrackHistos(s_trkFolder& folder, const Analysis::TauCommonDetails* p_taudetails, const Trk::TrackSummary* trackSummary, const Trk::Perigee* perigee);
  StatusCode fillCaloHistos(s_caloFolder& folder, Analysis::TauJetContainer::const_iterator ftau, const Analysis::TauCommonDetails* p_taudetails);
  StatusCode fillCombHistos(s_combFolder& folder,  const Analysis::TauCommonDetails* p_taudetails);
  StatusCode fillBDTHistos(s_BDTFolder&, Analysis::TauJetContainer::const_iterator ftau);
  StatusCode fillPhysicsHistograms(Analysis::TauJetContainer::const_iterator ftau, const Trk::Perigee* perigee);
  //StatusCode fillOldHistograms(Analysis::TauJetContainer::const_iterator ftau, const Analysis::TauCommonDetails* p_taudetails, const Analysis::TauPID* p_tauid);
  
  std::string fixName(std::string name);
  
};

#endif
