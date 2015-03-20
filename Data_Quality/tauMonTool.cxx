//*********************************************************************
// tauMonTool.cxx : Implementation of offline tau Data Quality
//                  histogramming, inheriting from ManagedMonitorToolBase
//
// authors:  C Cuenca Almenar, S Demers, E Ideal, A Leister, YALE
//*********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"

#include "tauMonitoring/tauMonTool.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#define CRACK_MIN 1.3
#define CRACK_MAX 1.7

#define PHIMIN -3.1415936
#define PHIMAX 3.1415936 

tauMonTool::tauMonTool( const std::string & type,
			const std::string & name, 
			const IInterface* parent ) :ManagedMonitorToolBase( type, name, parent ),
						    m_tauJetKey("TauRecContainer"),
						    m_tauCommonDetailsContainer("TauCommonDetailsContainer"),
						    m_doTrigger(false),
						    m_trigDec("Trig::TrigDecisionTool/TrigDecisionTool")
{
  declareProperty("doTrigger",
		  m_doTrigger,
		  "Run trigger aware monitoring. Only if trigger data is available");
  m_maxNLB = 1000;
}


tauMonTool::~tauMonTool() {}


//--------------------------------------------------------------------------------
// Book Histograms
//--------------------------------------------------------------------------------
StatusCode tauMonTool::bookHistograms() {
  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "StatusCode tauMonTool::bookHistograms, doTrigger = " << m_doTrigger << endreq;
    
    //--------------------
    // Top level folder
    //--------------------
    std::string folderName = "Tau";
    MonGroup mainFolder (this, folderName, run);
    
    if ( bookBasicPlots(m_basicPlots, mainFolder,"" ).isFailure() )
      log << MSG::ERROR << "Couldn't book basicPlots histograms" << endreq;

    
    //--------------------
    // Tau Barrel Monitoring groups      
    //--------------------
    if ( bookHistos(m_tauB, folderName + "/TauB", run) . isFailure())
      log << MSG::ERROR << "Couldn't book barrel histograms" << endreq;
    

    //--------------------
    // Tau Endcap Monitoring groups
    //--------------------
    if ( bookHistos(m_tauE, folderName + "/TauE", run) . isFailure())
      log << MSG::ERROR << "Couldn't book endcap histograms" << endreq;
    

    //--------------------
    // Tau Crack Monitoring groups
    //--------------------
    if ( bookHistos(m_tauCR, folderName + "/TauCR", run) . isFailure())
      log << MSG::ERROR << "Couldn't book crack histograms" << endreq;

    
    //--------------------
    // Tau Trigger-Aware Monitoring Groups
    //--------------------
    if (m_doTrigger) {
      if ( (m_trigDec.retrieve()).isSuccess()) {	
	
	//--------------------
	// create trigger folder and trigger summary histogram
	//--------------------
	MonGroup tau_Trigger (this, "Tau/Trigger", run);
	StatusCode sc = Book1DHist ( &m_triggers, &tau_Trigger, "m_triggers", "Fired triggers per tau candidate", 4, 0.5, 4.5);
	if (sc.isSuccess()) {
	  m_triggers->GetXaxis()->SetBinLabel(1,"EM");
	  m_triggers->GetXaxis()->SetBinLabel(2,"Jet");
	  m_triggers->GetXaxis()->SetBinLabel(3,"Taus");
	  m_triggers->GetXaxis()->SetBinLabel(4,"MBTS");
	} else {
	  log << MSG::ERROR << "Couldn't book m_triggers histogram" << endreq;
	}
	
	
	//--------------------
	// book trigger plots
	//--------------------
	MonGroup tau_Trigger_Jet (this, "Tau/Trigger/Jet" , run);
	MonGroup tau_Trigger_Em  (this, "Tau/Trigger/Em"  , run);
	MonGroup tau_Trigger_Tau (this, "Tau/Trigger/Tau" , run);
	MonGroup tau_Trigger_Mbts(this, "Tau/Trigger/Mbts", run);
	
	
	if ( bookBasicPlots(m_jetPlots,  tau_Trigger_Jet ,"jet_" ).isFailure() )
	  log << MSG::ERROR << "Couldn't book jetTrigger histograms" << endreq;
	
	if ( bookBasicPlots(m_emPlots,   tau_Trigger_Em  ,"em_"  ).isFailure() )
	  log << MSG::ERROR << "Couldn't book emTrigger histograms" << endreq;
	
	if ( bookBasicPlots(m_tauPlots,  tau_Trigger_Tau ,"tau_" ).isFailure() )
	  log << MSG::ERROR << "Couldn't book tauTrigger histograms" << endreq;
	
	if ( bookBasicPlots(m_mbtsPlots, tau_Trigger_Mbts,"mbts_").isFailure() )
	  log << MSG::ERROR << "Couldn't book mbtsTrigger histograms" << endreq;
	
      } else {
	log << MSG::ERROR << "couldn't retrieve trigger info" << endreq;
	return StatusCode::FAILURE;
      }    
    }
    
    //------------------------
    // book physics histograms
    //------------------------
    if (bookPhysicsHistograms() . isFailure() )
      log << MSG::ERROR << "Failed to book physics histograms" << endreq;
    
    //--------------------
    // book old histograms
    //--------------------
    /*if (bookOldHistograms() . isFailure() )
      log << MSG::ERROR << "Failed to book old histograms" << endreq;*/
  
  return StatusCode::SUCCESS;
}



//--------------------------------------------------------------------------------
// fill histograms
//--------------------------------------------------------------------------------
StatusCode tauMonTool::fillHistograms() {
  MsgStream log( msgSvc(), name() );
  //log << MSG::DEBUG << "StatusCode tauMonTool::fillHistograms" << endreq;
  log << MSG::DEBUG << "StatusCode tauMonTool::fillHistograms,  doTrigger = " << m_doTrigger << endreq;


  //--------------------
  //figure out current LB
  //--------------------
  const DataHandle<EventInfo> evtInfo;
  if ( (evtStore()->retrieve(evtInfo)) . isFailure()) {
    log << MSG::ERROR << "couldn't retrieve event info" << endreq;
    return StatusCode::FAILURE;
  }
 
  //LAr event veto: skip events rejected by LAr
  if(evtInfo->errorState(EventInfo::LAr)==EventInfo::Error){
    return StatusCode::SUCCESS;
  }

  EventID* eventID = evtInfo->event_ID();
  m_currentLB = eventID->lumi_block();
  
  

  //--------------------
  //Figure out trigger stuff
  //--------------------
  bool jetPassed  = false;
  bool emPassed   = false;
  bool tauPassed  = false;
  bool mbtsPassed = false;

  if (m_doTrigger) {
    if  (! m_trigDec.empty())  {
      const Trig::ChainGroup* jetChain  = m_trigDec->getChainGroup("L1_J.*");
      const Trig::ChainGroup* emChain   = m_trigDec->getChainGroup("L1_EM.*");
      const Trig::ChainGroup* tauChain  = m_trigDec->getChainGroup("L1_TAU.*");
      const Trig::ChainGroup* mbtsChain = m_trigDec->getChainGroup("L1_MBTS.*");
      
      jetPassed  = jetChain ->isPassed();
      emPassed   = emChain  ->isPassed();
      tauPassed  = tauChain ->isPassed();
      mbtsPassed = mbtsChain->isPassed();
      
      //how to get the list of items, ie:
      //std::vector<std::string>  jetList =  jetChain->getListOfTriggers();      
    } else {
      log << MSG::ERROR << "couldn't read trigger info" << endreq;
      return StatusCode::FAILURE;
    }    
  }
  
  //--------------------
  // Get StoreGate
  //--------------------
  if ( service("StoreGateSvc",m_storeGate).isFailure() ) {
    log << MSG::ERROR << "Unable to get pointer to StoreGateSvc" << endreq;
    return StatusCode::FAILURE;
  }



  //--------------------
  //Get Tau container
  //--------------------
  const Analysis::TauJetContainer  *tau_container;
  if ( (m_storeGate->retrieve(tau_container,m_tauJetKey)).isFailure() ) {
    log << MSG::WARNING << "Cannot retrieve " << m_tauJetKey << endreq;
    return StatusCode::FAILURE;
  } else {    

    int ntaus   = tau_container->size(); 
    int ntausB  = 0;
    int ntausCR = 0;
    int ntausE  = 0;
    int nHighPTtausB  = 0;
    int nHighPTtausCR = 0;
    int nHighPTtausE  = 0;


    //--------------------
    // Loop over taus (need beginning and end)
    //--------------------
    Analysis::TauJetContainer::const_iterator ftau = tau_container->begin();
    Analysis::TauJetContainer::const_iterator etau = tau_container->end();
    for (; ftau != etau; ftau++)  {
      
     
      //--------------------
      //prepare author vector
      //--------------------
      bool isCalo = (*ftau)->hasAuthor( TauJetParameters::tauRec);
      bool isTrk  = (*ftau)->hasAuthor( TauJetParameters::tau1P3P); 
      
      std::vector<int> author;
      if (isCalo)           author.push_back(1);
      if (isTrk)            author.push_back(2);
      if (isCalo && isTrk)  author.push_back(3);
      if (isCalo && !isTrk) author.push_back(4);
      if (!isCalo && isTrk) author.push_back(5);

      
      //--------------------
      // fill other basic histograms
      //--------------------
      fillBasicPlots( m_basicPlots,author,ftau);
      if (m_doTrigger) 
	if (! m_trigDec.empty()) {
	  if (  emPassed) { fillBasicPlots(m_emPlots,author,ftau);   m_triggers->Fill(1); }
	  if ( jetPassed) {	fillBasicPlots(m_jetPlots,author,ftau);  m_triggers->Fill(2); }
	  if ( tauPassed) {	fillBasicPlots(m_tauPlots,author,ftau);  m_triggers->Fill(3); }
	  if (mbtsPassed) {	fillBasicPlots(m_mbtsPlots,author,ftau); m_triggers->Fill(4); }
	}
      
      //--------------------
      // ID info
      //--------------------
      const Analysis::TauPID* p_tauid = (*ftau)->tauID();

      //--------------------
      // track info
      //--------------------
      const Trk::TrackSummary *trackSummary = NULL;
      const Trk::Perigee *perigee = NULL;

      double pTGev = 0.0;
      
      if ((*ftau)->numTrack() >= 1) {
	trackSummary = (*ftau)->track(0)->trackSummary();
	perigee = (*ftau)->track(0)->perigee();
	pTGev = perigee->pT()/1000.0;
      }


      //--------------------
      // fill histograms
      //--------------------

      //get details
      const Analysis::TauCommonDetails* p_taudetails = (*ftau)->details<const Analysis::TauCommonDetails>();

      float abs_eta = fabs((*ftau)->eta());

      if (abs_eta <= CRACK_MIN ) {
	if (fillHistograms(m_tauB, author, ftau, p_taudetails, p_tauid, trackSummary, perigee).isFailure()) log<<MSG::ERROR<<"Failed to fill barrel histograms" 
													       <<endreq;
	ntausB++;
	if (pTGev > 100.0){ nHighPTtausB++; }
      }	  
      if (abs_eta > CRACK_MAX ) {
	if (fillHistograms(m_tauE, author, ftau, p_taudetails, p_tauid, trackSummary, perigee).isFailure()) log<<MSG::ERROR<<"Failed to fill endcap histograms" 
													       <<endreq;
	ntausE++;
	if (pTGev > 100.0){ nHighPTtausE++; }
      }
      if (abs_eta > CRACK_MIN && abs_eta < CRACK_MAX ) {
	if (fillHistograms(m_tauCR, author, ftau, p_taudetails, p_tauid, trackSummary, perigee).isFailure()) log<<MSG::ERROR<<"Failed to fill crack histograms" 
														<<endreq;
	ntausCR++;
	if (pTGev > 100.0){ nHighPTtausCR++; }
      }

      
      //-----------------------
      //Fill physics histograms
      //-----------------------
      if(fillPhysicsHistograms(ftau, perigee) . isFailure())
	log << MSG::ERROR << "Failed to fill physics histograms" << endreq;


      //--------------------
      // fill old histograms
      //--------------------
      /*if (fillOldHistograms(ftau, p_taudetails, p_tauid) . isFailure())
	log << MSG::ERROR << "Failed to fill old histograms" << endreq;*/
      

    } // end of loop over taus
    

    //--------------------
    // fill ntaus
    //--------------------
    m_basicPlots.h_ntaus->Fill(ntaus);
    m_basicPlots.h_nHighPTtaus->Fill(nHighPTtausB + nHighPTtausCR + nHighPTtausE);
    if (m_doTrigger)
      if (! m_trigDec.empty()) {
	if ( jetPassed) m_jetPlots .h_ntaus->Fill(ntaus);
	if (  emPassed) m_emPlots  .h_ntaus->Fill(ntaus);
	if ( tauPassed) m_tauPlots .h_ntaus->Fill(ntaus);
	if (mbtsPassed) m_mbtsPlots.h_ntaus->Fill(ntaus);
      }
    
    m_tauB .kinFolder.h_ntaus->Fill(ntausB);
    m_tauCR.kinFolder.h_ntaus->Fill(ntausCR);
    m_tauE .kinFolder.h_ntaus->Fill(ntausE);
    m_tauB .trkFolder.h_nHighPTtaus->Fill(nHighPTtausB);
    m_tauCR.trkFolder.h_nHighPTtaus->Fill(nHighPTtausCR);
    m_tauE .trkFolder.h_nHighPTtaus->Fill(nHighPTtausE);

  }//if ( (m_storeGate->retrieve(tau_container
  
  return StatusCode::SUCCESS;  
}


//--------------------------------------------------------------------------------
// book and register a 1D histogram
//--------------------------------------------------------------------------------
StatusCode tauMonTool::Book1DHist (TH1 ** hist, MonGroup * monName, std::string  histName, std::string  histTitle, int NBins, double lowBin, double highBin)
{

  *hist = new TH1F( histName.c_str(), histTitle.c_str(), NBins, lowBin, highBin);
  StatusCode sc = monName->regHist( *hist);
  if ( sc.isFailure() )
    {
      MsgStream log( msgSvc(), name() );
      log << MSG::FATAL << "Failed to register " << histName << endreq;
      return sc;
    }
  return sc;
}


//--------------------------------------------------------------------------------
// book and register a 2D histogram
//--------------------------------------------------------------------------------
StatusCode tauMonTool::Book2DHist (TH2 ** hist, MonGroup * monName, std::string  histName, std::string  histTitle, int NXBins, double lowXBin, double highXBin, int NYBins, double lowYBin, double highYBin)
{

  *hist = new TH2F( histName.c_str(), histTitle.c_str(), NXBins, lowXBin, highXBin, NYBins, lowYBin, highYBin);
  StatusCode sc = monName->regHist( *hist);
  if ( sc.isFailure() )
    {
      MsgStream log( msgSvc(), name() );
      log << MSG::FATAL << "Failed to register " << histName << endreq;
      return sc;
    }
  return sc;
}
 

//--------------------------------------------------------------------------------
// post processing
//--------------------------------------------------------------------------------
StatusCode tauMonTool::procHistograms()
{

  if( endOfRun || endOfLumiBlock) {

    //check that we have enough bins
    //if  (m_currentLB > m_maxNLB) {
      //m_maxNLB = 2*m_maxNLB;
      //for (int i=0; i<6; i++) m_ntaus_vs_LB[i]->GetXaxis()->SetRangeUser(0,m_maxNLB);
      
    
    //fill the vs_LB plots
    for (unsigned int iN=0; iN<6; iN++) {
      if (m_currentLB <= m_maxNLB) {
	m_basicPlots.h_ntaus_vs_LB[iN]->SetBinContent(m_currentLB,m_basicPlots.h_ntausLB[iN]);
	if (m_doTrigger) 
	  if  (! m_trigDec.empty()) {
	    m_jetPlots  .h_ntaus_vs_LB[iN]->SetBinContent(m_currentLB,  m_jetPlots.h_ntausLB[iN]);
	    m_emPlots   .h_ntaus_vs_LB[iN]->SetBinContent(m_currentLB,   m_emPlots.h_ntausLB[iN]);
	    m_tauPlots  .h_ntaus_vs_LB[iN]->SetBinContent(m_currentLB,  m_tauPlots.h_ntausLB[iN]);
	    m_mbtsPlots .h_ntaus_vs_LB[iN]->SetBinContent(m_currentLB, m_mbtsPlots.h_ntausLB[iN]);
	  }
      }

      //reset counters
      m_basicPlots.h_ntausLB[iN] = 0;
      if (m_doTrigger) 
	if  (! m_trigDec.empty()) {
	  m_jetPlots  .h_ntausLB[iN] = 0;
	  m_emPlots   .h_ntausLB[iN] = 0;
	  m_tauPlots  .h_ntausLB[iN] = 0;
	  m_mbtsPlots .h_ntausLB[iN] = 0;    
	}
    }
  }
  
  return StatusCode::SUCCESS;
  
}

 
StatusCode tauMonTool::bookBasicPlots(s_basicPlots& someBasicPlots, MonGroup &aGroup, std::string prefix){
  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "StatusCode tauMonTool::bookBasicPlots, " << prefix << endreq;


  for (unsigned int iN=0; iN<6; iN++)
    someBasicPlots.h_ntausLB[iN]=0;

  std::string name = prefix + "tauAuthor";
  StatusCode sc = Book1DHist ( &someBasicPlots.h_author, &aGroup, name, "Author of tau candidates", 5, 0.5, 5.5 );

  someBasicPlots.h_author->GetXaxis()->SetBinLabel(1, "All Caloseed");
  someBasicPlots.h_author->GetXaxis()->SetBinLabel(2, "All Trackseed");
  someBasicPlots.h_author->GetXaxis()->SetBinLabel(3, "Both Seeds");
  someBasicPlots.h_author->GetXaxis()->SetBinLabel(4, "Caloseed ONLY");
  someBasicPlots.h_author->GetXaxis()->SetBinLabel(5, "Trackseed ONLY");

  name = prefix + "nTauCandidates";
  sc = Book1DHist ( &someBasicPlots.h_ntaus, &aGroup, name, "Number of tau candidates", 31, -0.5, 30.5);
  someBasicPlots.h_ntaus->GetXaxis()->SetTitle("Number of Taus per Event");

  name = prefix + "tauEta";
  sc = Book1DHist ( &someBasicPlots.h_eta, &aGroup, name, "Eta of tau candidates", 51, -2.55, 2.55);
  someBasicPlots.h_eta->GetXaxis()->SetTitle("Eta");
  someBasicPlots.h_eta->GetYaxis()->SetTitle("Number of Candidates");

  name = prefix + "tauPhi";
  sc = Book1DHist ( &someBasicPlots.h_phi, &aGroup, name, "Phi of tau candidates", 65, -3.1415936-0.098174/2., 3.1415936+0.098174/2.);
  someBasicPlots.h_phi->GetXaxis()->SetTitle("Phi");
  someBasicPlots.h_phi->GetYaxis()->SetTitle("Number of Candidates");

  name = prefix + "tauEt";
  sc = Book1DHist ( &someBasicPlots.h_et, &aGroup, name, "Et of tau candidates", 50, 0.0, 250.0);
  someBasicPlots.h_et->GetXaxis()->SetTitle("Transverse Energy (GeV)");
  someBasicPlots.h_et->GetYaxis()->SetTitle("Number of Candidates");

  name = prefix + "tauCharge";
  sc = Book1DHist ( &someBasicPlots.h_charge, &aGroup, name, "Charge of tau candidates", 11, -5.5, 5.5);
  someBasicPlots.h_charge->GetXaxis()->SetTitle("Charge");
  someBasicPlots.h_charge->GetYaxis()->SetTitle("Number of Candidates");

  name = prefix + "tauNumTracks";
  sc = Book1DHist ( &someBasicPlots.h_numTracks, &aGroup, name, "Number of Tracks for tau candidates", 21, -0.5, 20.5);
  someBasicPlots.h_numTracks->GetXaxis()->SetTitle("Number of Tracks");
  someBasicPlots.h_numTracks->GetYaxis()->SetTitle("Number of Candidates");

  name = prefix + "nHighPtTauCandidates";
  sc = Book1DHist ( &someBasicPlots.h_nHighPTtaus, &aGroup, name, "Number of High pT tau candidates", 31, -0.5, 30.5);
  someBasicPlots.h_ntaus->GetXaxis()->SetTitle("Number of Taus per Event");
  
  name = prefix + "tauEtVsEta";
  sc = Book2DHist ( &someBasicPlots.h_EtVsEta, &aGroup, name, "Tau Et vs. Eta",  51, -2.55, 2.55, 100, 0, 200);
  someBasicPlots.h_EtVsEta->GetXaxis()->SetTitle("Eta");
  someBasicPlots.h_EtVsEta->GetYaxis()->SetTitle("Transverse Energy (GeV)");
  someBasicPlots.h_EtVsEta->GetZaxis()->SetTitle("Number of Candidates");


  name = prefix + "tauEtVsPhi";
  sc = Book2DHist ( &someBasicPlots.h_EtVsPhi, &aGroup, name, "Tau Et vs. Phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 100, 0, 200);
  someBasicPlots.h_EtVsPhi->GetXaxis()->SetTitle("Phi");
  someBasicPlots.h_EtVsPhi->GetYaxis()->SetTitle("Transverse Energy (GeV)");


  name = prefix + "tauPhiVsEta";
  sc = Book2DHist ( &someBasicPlots.h_PhiVsEta, &aGroup, name, "Tau Phi vs. Eta", 51, -2.55, 2.55, 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64.);
  someBasicPlots.h_PhiVsEta->GetXaxis()->SetTitle("Eta");
  someBasicPlots.h_PhiVsEta->GetYaxis()->SetTitle("Phi");
  someBasicPlots.h_PhiVsEta->GetZaxis()->SetTitle("Number of Candidates");
  

  name = prefix + "tauEtaVsLB";
  sc = Book2DHist ( &someBasicPlots.h_Eta_vs_LB, &aGroup, name, "Tau Eta vs Lumiblock", 51, -2.55, 2.55, m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0 );
  someBasicPlots.h_Eta_vs_LB->GetXaxis()->SetTitle("Eta");
  someBasicPlots.h_Eta_vs_LB->GetYaxis()->SetTitle("Lumiblock");


  name = prefix + "tauPhiVsLB";
  sc = Book2DHist ( &someBasicPlots.h_Phi_vs_LB, &aGroup, name, "Tau Phi vs Lumiblock", 65, PHIMIN+PHIMIN/64, PHIMAX+PHIMAX/64, m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0 );
  someBasicPlots.h_Phi_vs_LB->GetXaxis()->SetTitle("Phi");
  someBasicPlots.h_Phi_vs_LB->GetYaxis()->SetTitle("Lumiblock");


  name = prefix + "z_ntauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[0], &aGroup, name, "Total number of tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[0]->GetXaxis()->SetTitle("Luminosity Block");
  someBasicPlots.h_ntaus_vs_LB[0]->GetYaxis()->SetTitle("Number of Candidates");


  name = prefix + "zz_nCaloTauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[1], &aGroup, name, "Number of calo-seeded tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[1]->GetXaxis()->SetTitle("Luminosity Block");
  someBasicPlots.h_ntaus_vs_LB[1]->GetYaxis()->SetTitle("Number of Candidates");


  name = prefix + "zz_nTrkTauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[2], &aGroup, name, "Number of trk-seeded tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[2]->GetXaxis()->SetTitle("Luminosity Block");


  name = prefix + "zzz_nMergTauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[3], &aGroup, name, "Number of merged tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[3]->GetXaxis()->SetTitle("Luminosity Block");
  someBasicPlots.h_ntaus_vs_LB[3]->GetYaxis()->SetTitle("Number of Candidates");


  name = prefix + "zzzz_nXCaloTauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[4], &aGroup, name, "Number of excl cal tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[4]->GetXaxis()->SetTitle("Luminosity Block");


  name = prefix + "zzzz_nXTrkTauLB";
  sc = Book1DHist ( &someBasicPlots.h_ntaus_vs_LB[5], &aGroup, name, "Number of excl trk tau candidates per LB", m_maxNLB, 0, m_maxNLB);
  someBasicPlots.h_ntaus_vs_LB[5]->GetXaxis()->SetTitle("Luminosity Block");
  someBasicPlots.h_ntaus_vs_LB[5]->GetYaxis()->SetTitle("Number of Candidates");
  
  return StatusCode::SUCCESS;
}


 void tauMonTool::fillBasicPlots(s_basicPlots& someBasicPlots, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau){

  someBasicPlots.h_ntausLB[0]++;
  
  for (unsigned int iA=0; iA < author.size(); iA++) {
    someBasicPlots.h_author->Fill(author[iA]);
    someBasicPlots.h_ntausLB[author[iA]]++;
  }


  float eta     = (*ftau)->eta();
  float et      = (*ftau)->et() / 1000.0;
  float phi     = (*ftau)->phi();
  int   charge  = (int) (*ftau)->charge();
  int numTracks = (int) (*ftau)->numTrack();

  someBasicPlots.h_eta      ->Fill( eta );
  someBasicPlots.h_phi      ->Fill( phi );
  someBasicPlots.h_et       ->Fill( et );
  someBasicPlots.h_charge   ->Fill( charge);
  someBasicPlots.h_numTracks->Fill( numTracks );

  someBasicPlots.h_EtVsEta  ->Fill( eta, et);
  someBasicPlots.h_EtVsPhi  ->Fill( phi, et);
  someBasicPlots.h_PhiVsEta ->Fill( eta, phi);

  someBasicPlots.h_Eta_vs_LB->Fill( eta, m_currentLB);
  someBasicPlots.h_Phi_vs_LB->Fill( phi, m_currentLB);
  
  return;
}


StatusCode tauMonTool::bookHistos(s_mainFolder& mainFolder, std::string folderName, Interval_t interval){
  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "StatusCode tauMonTool::bookHistos, folderName = " << folderName << endreq;

  MonGroup folder(this, folderName, interval);
  
  if ( bookKinHistos(mainFolder.kinFolder, folder) .isFailure() )
    log << MSG::ERROR << "Couldn't book kinematic histograms" << endreq;
  
  if ( bookIDHistos(mainFolder.idFolder, folderName, interval).isFailure() )
    log << MSG::ERROR << "Couldn't book identification histograms" << endreq;
  
  if ( bookTrackHistos(mainFolder.trkFolder, folderName, interval).isFailure() )
    log << MSG::ERROR << "Couldn't book track histograms" << endreq;
  
  if ( bookCaloHistos(mainFolder.caloFolder, folderName, interval).isFailure() )
    log << MSG::ERROR << "Couldn't book calorimeter histograms" << endreq;
  
  if ( bookCombHistos(mainFolder.combFolder, folderName, interval).isFailure() )
    log << MSG::ERROR << "Couldn't book combined histograms" << endreq;
  
  return StatusCode::SUCCESS;
}
  
StatusCode tauMonTool::bookKinHistos(s_kinFolder& folder,  MonGroup &aGroup){

  StatusCode sc = StatusCode::SUCCESS; 

  std::string prefix = this->fixName(aGroup.system());  

  std::string name = prefix + "_tauAuthor";
  sc = Book1DHist ( &folder.h_author, &aGroup, name, "Author of tau candidate", 5, 0.5, 5.5 );

  folder.h_author->GetXaxis()->SetBinLabel(1, "All Caloseed");
  folder.h_author->GetXaxis()->SetBinLabel(2, "All Trackseed");
  folder.h_author->GetXaxis()->SetBinLabel(3, "Both Seeds");
  folder.h_author->GetXaxis()->SetBinLabel(4, "Caloseed ONLY");
  folder.h_author->GetXaxis()->SetBinLabel(5, "Trackseed ONLY");

  name = prefix + "_nTauCandidates";
  sc = Book1DHist ( &folder.h_ntaus, &aGroup, name, "Number of tau candidates", 31, -0.5, 30.5);
  folder.h_ntaus->GetXaxis()->SetTitle("Number of Taus per Event");

  name = prefix + "_tauEta";
  sc = Book1DHist ( &folder.h_eta, &aGroup, name, "Eta of tau candidates", 51, -2.55, 2.55);
  folder.h_eta->GetXaxis()->SetTitle("Eta");
  folder.h_eta->GetYaxis()->SetTitle("Number of Candidates");


  name = prefix + "_tauPhi";
  sc = Book1DHist ( &folder.h_phi, &aGroup, name, "Phi of tau candidates", 65, -3.1415936-0.098174/2., 3.1415936+0.098174/2.);
  folder.h_phi->GetXaxis()->SetTitle("Phi");


  name = prefix + "_tauEt";
  sc = Book1DHist ( &folder.h_et, &aGroup, name, "Et of tau candidates", 50, 0.0, 250.0);
  folder.h_et->GetXaxis()->SetTitle("Transverse Energy (GeV)");


  name = prefix + "_tauCharge";
  sc = Book1DHist ( &folder.h_charge, &aGroup, name, "Charge of tau candidates", 11, -5.5, 5.5);
  folder.h_charge->GetXaxis()->SetTitle("Charge");
  folder.h_charge->GetYaxis()->SetTitle("Number of Candidates");


  name = prefix + "_tauPhiVsEta";
  sc = Book2DHist ( &folder.h_PhiVsEta, &aGroup, name, "Tau Phi vs. Eta", 51, -2.55, 2.55, 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64.);
  folder.h_PhiVsEta->GetXaxis()->SetTitle("Eta");
  folder.h_PhiVsEta->GetYaxis()->SetTitle("Phi");
  folder.h_PhiVsEta->GetZaxis()->SetTitle("Number of Candidates");
  

  name = prefix + "_tauEtaVsLB";
  sc = Book2DHist ( &folder.h_Eta_vs_LB, &aGroup, name, "Tau Eta vs Lumiblock", 51, -2.55, 2.55, m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0);
  folder.h_Eta_vs_LB->GetXaxis()->SetTitle("Eta");
  folder.h_Eta_vs_LB->GetYaxis()->SetTitle("Lumiblock");


  name = prefix + "_tauPhiVsLB";
  sc = Book2DHist ( &folder.h_Phi_vs_LB, &aGroup, name, "Tau Phi vs Lumiblock", 65, PHIMIN+PHIMIN/64, PHIMAX+PHIMAX/64, m_maxNLB/10+1, -0.5, (double)m_maxNLB+0.5);
  folder.h_Phi_vs_LB->GetXaxis()->SetTitle("Phi");
  folder.h_Phi_vs_LB->GetYaxis()->SetTitle("Lumiblock");

  return sc;
}


std::string tauMonTool::fixName(std::string name) {
  
  std::string::size_type start = 0;
  
  while ( (start = name.find('/')) != std::string::npos)
    name.replace(start,1,"_");
  
  return name;
}
  
StatusCode tauMonTool::bookIDHistos(s_idFolder& folder,std::string folderName, Interval_t interval) {

  folderName = folderName + "/Identification";
  MonGroup aGroup(this, folderName, interval);

  if ( bookBDTLooseHistos(folder.BDTLooseFolder, folderName, interval).isFailure() ){
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Couldn't book BDTLoose histograms" << endreq;
  }

  if ( bookBDTMedHistos(folder.BDTMedFolder, folderName, interval).isFailure() ){
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Couldn't book BDTMedium histograms" << endreq;
  }

  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;

  sc = Book1DHist ( &folder.h_tauCutLoose,          &aGroup, folderName + "_tauCutLoose",          "Identification Flag: tauCutLoose",          2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauCutLoose->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauCutLoose->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_tauCutMedium,         &aGroup, folderName + "_tauCutMedium",         "Identification Flag: tauCutMedium",         2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauCutMedium->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauCutMedium->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_tauCutTight,          &aGroup, folderName + "_tauCutTight",          "Identification Flag: tauCutTight",          2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauCutTight->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauCutTight->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_tauLlhLoose,          &aGroup, folderName + "_tauLlhLoose",          "Identification Flag: tauLlhLoose",          2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauLlhLoose->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauLlhLoose->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_tauLlhMedium,         &aGroup, folderName + "_tauLlhMedium",         "Identification Flag: tauLlhMedium",         2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauLlhMedium->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauLlhMedium->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_tauLlhTight,          &aGroup, folderName + "_tauLlhTight",          "Identification Flag: tauLlhTight",          2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_tauLlhTight->GetXaxis()->SetBinLabel(1, "False");
    folder.h_tauLlhTight->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_electronVetoLoose,    &aGroup, folderName + "_electronVetoLoose",    "Loose electron Veto",                       2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_electronVetoLoose->GetXaxis()->SetBinLabel(1, "False");
    folder.h_electronVetoLoose->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_electronVetoMedium,   &aGroup, folderName + "_electronVetoMedium",   "Medium electron Veto",                      2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_electronVetoMedium->GetXaxis()->SetBinLabel(1, "False");
    folder.h_electronVetoMedium->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_electronVetoTight,    &aGroup, folderName + "_electronVetoTight",    "Tight electron Veto",                       2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_electronVetoTight->GetXaxis()->SetBinLabel(1, "False");
    folder.h_electronVetoTight->GetXaxis()->SetBinLabel(2, "True");  
  }
  
  sc = Book1DHist ( &folder.h_muonVeto,             &aGroup, folderName + "_muonVeto",             "Muon Veto",                                 2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_muonVeto->GetXaxis()->SetBinLabel(1, "False");
    folder.h_muonVeto->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_eleBDTLoose,          &aGroup, folderName + "_eleBDTLoose",          "Loose BDT",                                 2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_eleBDTLoose->GetXaxis()->SetBinLabel(1, "False");
    folder.h_eleBDTLoose->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_eleBDTMedium,         &aGroup, folderName + "_eleBDTMedium",         "Medium BDT",                                2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_eleBDTMedium->GetXaxis()->SetBinLabel(1, "False");
    folder.h_eleBDTMedium->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_eleBDTTight,          &aGroup, folderName + "_eleBDTTight",          "Tight BDT",                                 2, -0.5, 1.5);
  if (sc.isSuccess()) {
    folder.h_eleBDTTight->GetXaxis()->SetBinLabel(1, "False");
    folder.h_eleBDTTight->GetXaxis()->SetBinLabel(2, "True");  
  }

  sc = Book1DHist ( &folder.h_BDTJetScore,         &aGroup, folderName + "_BDTJetScore",         "BDT Score for Jet Rejection",                             48, -0.1, 1.1);
  if (sc.isSuccess()) 
      folder.h_BDTJetScore->GetXaxis()->SetTitle("Boosted Decision Tree Score"); 
  
  sc = Book1DHist ( &folder.h_BDTEleScore,         &aGroup, folderName + "_BDTEleScore",         "BDT Score for Electron Rejection",                        48, -0.1, 1.1);
  if (sc.isSuccess()) 
      folder.h_BDTEleScore->GetXaxis()->SetTitle("Boosted Decision Tree Score");
  
  sc = Book1DHist ( &folder.h_BDTJetScoreSigTrans, &aGroup, folderName + "_BDTJetScoreSigTrans", "Signal Transformed BDT Score for Jet Rejection ",         48, -0.1, 1.1);
  if (sc.isSuccess()) 
      folder.h_BDTJetScoreSigTrans->GetXaxis()->SetTitle("Boosted Decision Tree Score");

  sc = Book1DHist ( &folder.h_BDTJetScoreBkgTrans, &aGroup, folderName + "_BDTJetScoreBkgTrans", "Background Transformed BDT Score for Electron Rejection", 48, -0.1, 1.1);
  if (sc.isSuccess()) 
      folder.h_BDTJetScoreBkgTrans->GetXaxis()->SetTitle("Boosted Decision Tree Score");

  return sc;
}

 StatusCode tauMonTool::bookTrackHistos(s_trkFolder& folder,std::string folderName, Interval_t interval) {


  folderName = folderName + "/Track";
  MonGroup aGroup(this, folderName, interval);
  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;
  
  sc = Book1DHist ( &folder.h_leadTrkPt, &aGroup, folderName + "_leadTrkPt","Pt of Leading track", 50, 0, 150);
  if (sc.isSuccess()) folder.h_leadTrkPt->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  
  sc = Book1DHist ( &folder.h_massTrkSys, &aGroup, folderName + "_massTrkSys","Mass of the Track System", 30, -1.0, 5.0);
  if (sc.isSuccess()) folder.h_massTrkSys->GetXaxis()->SetTitle("Invariant Mass (GeV)");

  sc = Book1DHist ( &folder.h_trkWidth2, &aGroup, folderName + "_trkWidth2","Weighted Track Width", 25, 0.0, 0.1);
  if (sc.isSuccess()) folder.h_trkWidth2->GetXaxis()->SetTitle("Momentum-Weighted Width of Track System");

  sc = Book1DHist ( &folder.h_trFlightPathSig, &aGroup, folderName + "_trFlightPathSig","Track Transverse Flight Path Significance", 50, -5.0, 5.0);
  if (sc.isSuccess()) folder.h_trFlightPathSig->GetXaxis()->SetTitle("Transverse Flight Path Significance");

  sc = Book1DHist ( &folder.h_ipSigLeadTrk, &aGroup, folderName + "_ipSigLeadTrk","Impact Parameter Significance of Leading Track", 50, -5.0, 5.0);
  if (sc.isSuccess()) folder.h_ipSigLeadTrk->GetXaxis()->SetTitle("Transverse Impact Parameter Significance");

  sc = Book1DHist ( &folder.h_ipZ0SinThetaSigLeadTrk, &aGroup, folderName + "_ipZ0SinThetaSigLeadTrk","Impact Parameter z0 Sine Theta Significance of Leading Track", 50, -10.0, 10.0);
  if (sc.isSuccess()) folder.h_ipZ0SinThetaSigLeadTrk->GetXaxis()->SetTitle("Z0SinTheta Significance");

  sc = Book1DHist ( &folder.h_d0, &aGroup, folderName + "_d0","Track d0", 50, -5.0, 15.0);
  if (sc.isSuccess()) folder.h_d0->GetXaxis()->SetTitle("Transverse Impact Parameter (mm)");

  sc = Book1DHist ( &folder.h_z0, &aGroup, folderName + "_z0","Track z0", 50, -25.0, 25.0);
  if (sc.isSuccess()) folder.h_z0->GetXaxis()->SetTitle("Longitudinal Impact Parameter (mm)");

  sc = Book1DHist ( &folder.h_phi, &aGroup, folderName + "_phi","Track Phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64.);
  if (sc.isSuccess()) folder.h_phi->GetXaxis()->SetTitle("Phi");

  sc = Book1DHist ( &folder.h_eta, &aGroup, folderName + "_eta","Track Eta", 51, -2.55, 2.55);
  if (sc.isSuccess()) folder.h_eta->GetXaxis()->SetTitle("Eta");

  sc = Book1DHist ( &folder.h_pT, &aGroup, folderName + "_pT","Track pT", 50, 0.0, 150.0);
  if (sc.isSuccess()) folder.h_pT->GetXaxis()->SetTitle("Transverse Momentum (GeV)");

  sc = Book1DHist ( &folder.h_nHighPTtaus, &aGroup, folderName + "_nHighPTtaus","Number of High-pT Tau Candidates", 31, -0.5, 30.5);
  if (sc.isSuccess()) folder.h_nHighPTtaus->GetXaxis()->SetTitle("Number of Taus per Event");

  sc = Book1DHist ( &folder.h_numberOfTRTHighThresholdHits, &aGroup, folderName + "_numberOfTRTHighThresholdHits","Number of TRT High Threshold Hits", 15, -0.5, 14.5);
  if (sc.isSuccess()) folder.h_numberOfTRTHighThresholdHits->GetXaxis()->SetTitle("Number of High Threshold TRT Hits");

  sc = Book1DHist ( &folder.h_numberOfTRTHits, &aGroup, folderName + "_numberOfTRTHits","Number of TRT Low Threshold Hits", 101, -0.5, 100.5);
  if (sc.isSuccess()) folder.h_numberOfTRTHits->GetXaxis()->SetTitle("Number of Low Threshold TRT Hits");

  sc = Book1DHist ( &folder.h_numberOfTRTHighThresholdOutliers, &aGroup, folderName + "_numberOfTRTHighThresholdOutliers","Number of TRT High Threshold Outliers", 26, -0.5, 25.5);
  if (sc.isSuccess()) folder.h_numberOfTRTHighThresholdOutliers->GetXaxis()->SetTitle("Number of TRT High Threshold Outliers");

  sc = Book1DHist ( &folder.h_numberOfTRTOutliers, &aGroup, folderName + "_numberOfTRTOutliers","Number of TRT Low Threshold Outliers", 26, -0.5, 25.5);
  if (sc.isSuccess()) folder.h_numberOfTRTOutliers->GetXaxis()->SetTitle("Number of TRT Low Threshold Outliers");

  sc = Book1DHist ( &folder.h_numberOfSCTHits, &aGroup, folderName + "_numberOfSCTHits","Number of SCT Hits", 51, -0.5, 50.5);
  if (sc.isSuccess()) folder.h_numberOfSCTHits->GetXaxis()->SetTitle("Number of SCT Hits");

  sc = Book1DHist ( &folder.h_numberOfPixelHits, &aGroup, folderName + "_numberOfPixelHits","Number of Pixel Hits", 26, -0.5, 25.5);
  if (sc.isSuccess()) folder.h_numberOfPixelHits->GetXaxis()->SetTitle("Number of Pixel Hits");

  
  sc= Book2DHist (&folder.h_z0_vs_LB, &aGroup, folderName + "_z0VsLB", "Track z0 vs Lumiblock" , 50, -25.0, 25.0 , m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0);
  if (sc.isSuccess()) {
    folder.h_z0_vs_LB->GetXaxis()->SetTitle("Longitudinal Impact Parameter (mm)");
    folder.h_z0_vs_LB->GetYaxis()->SetTitle("Lumiblock");
  }

  return sc;
}


 StatusCode tauMonTool::bookCaloHistos(s_caloFolder& folder,std::string folderName, Interval_t interval) {


  folderName = folderName + "/Calo";
  MonGroup aGroup(this, folderName, interval);
  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;
  

  sc = Book1DHist ( &folder.h_eta, &aGroup, folderName + "_eta","Calorimeter eta of tau candidates", 51, -2.55, 2.55);
  if (sc.isSuccess()) {
    folder.h_eta->GetXaxis()->SetTitle("Eta");
    folder.h_eta->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_phi, &aGroup, folderName + "_phi","Calorimeter phi of tau candidates", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64.);
  if (sc.isSuccess()) {
    folder.h_phi->GetXaxis()->SetTitle("Phi");
    folder.h_phi->GetYaxis()->SetTitle("Number of Candidates");
  
  }

  sc = Book1DHist ( &folder.h_etEMCalib, &aGroup, folderName + "_etEMCalib","Calibrated EM ET of tau candidates", 50, 0.0, 150.0);
  if (sc.isSuccess()) {
    folder.h_etEMCalib->GetXaxis()->SetTitle("Calibrated EM ET (in GeV)");
    folder.h_etEMCalib->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_etHadCalib, &aGroup, folderName + "_etHadCalib","Calibrated Had ET of tau candidates", 50, 0.0, 150.0);
  if (sc.isSuccess()) {
    folder.h_etHadCalib->GetXaxis()->SetTitle("Calibrated Hadronic ET (in GeV)");
    folder.h_etHadCalib->GetYaxis()->SetTitle("Number of Candidates");
  }

  
  sc = Book1DHist ( &folder.h_etHadAtEMScale, &aGroup, folderName + "_etHadAtEMScale","Hadronic Energy at the EM Scale", 50, 0.0, 150.0);
  if (sc.isSuccess()) {
    folder.h_etHadAtEMScale->GetXaxis()->SetTitle("Had Et (GeV)");
    folder.h_etHadAtEMScale->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_etEMAtEMScale, &aGroup, folderName + "_etEMAtEMScale","EM energy at the EM scale", 50, 0.0, 150.0);
  if (sc.isSuccess()) {
    folder.h_etEMAtEMScale->GetXaxis()->SetTitle("EM Et (GeV)");
    folder.h_etEMAtEMScale->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_EMRadius, &aGroup, folderName + "_EMRadius","Uncalibrated EM Radius", 50, 0.0, 1.0);
  if (sc.isSuccess()) {
    folder.h_EMRadius->GetXaxis()->SetTitle("EM Radius");
    folder.h_EMRadius->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_centFrac, &aGroup, folderName + "_centralityFraction","Centrality fraction of tau candidates", 51, 0.0, 1.02);
  if (sc.isSuccess()) {
    folder.h_centFrac->GetXaxis()->SetTitle("Centrality Fraction");
    folder.h_centFrac->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_hadRadius, &aGroup, folderName + "_hadRadius","Hadronic Radius of tau candidates", 50, 0.0, 1.0);
  if (sc.isSuccess()) {    
    folder.h_hadRadius->GetXaxis()->SetTitle("Hadronic Radius");
    folder.h_hadRadius->GetYaxis()->SetTitle("Number of Candidiates");
  }

  sc = Book1DHist ( &folder.h_isolFrac, &aGroup, folderName + "_isolFrac","Isolation Fraction", 51, 0.0, 1.02);
  if (sc.isSuccess()) {    
    folder.h_isolFrac->GetXaxis()->SetTitle("Et Isolation Fraction");
    folder.h_isolFrac->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_stripWidth2, &aGroup, folderName + "_stripWidth2","Strip Width of tau candidates", 50, -0.1, 0.1);
  if (sc.isSuccess()) {    
    folder.h_stripWidth2->GetXaxis()->SetTitle("Strip Width");
    folder.h_stripWidth2->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_nStrip, &aGroup, folderName + "_nStrip","Number of strip cells of tau candidates", 50, -0.5, 49.5);
  if (sc.isSuccess()) {    
    folder.h_nStrip->GetXaxis()->SetTitle("Number of Strip Cells");
    folder.h_nStrip->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_trkAvgDist, &aGroup, folderName + "_trkAvgDist","Average Track Distance from Calorimeter Seed", 50, -1.0, 1.0);
  if (sc.isSuccess()){
    folder.h_trkAvgDist->GetXaxis()->SetTitle("Distance (mm)");
    folder.h_trkAvgDist->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_dRmax, &aGroup, folderName + "_dRmax", "Maximum delta R", 42, 0.0, 0.21);
  if (sc.isSuccess()){
    folder.h_dRmax->GetXaxis()->SetTitle("delta R");
    folder.h_dRmax->GetYaxis()->SetTitle("Number of Candidates");
  }

  /*
  sc = Book1DHist ( &folder.h_cellBasedEnergyRing[0], &aGroup, folderName + "_cellBasedEnergyRing1", "Cell Based Energy in Ring 1: 0 < R < 0.05", 50, 0, 100);
  if (sc.isSuccess()){
    folder.h_dRmax->GetXaxis()->SetTitle("Energy (GeV)");
    folder.h_dRmax->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_cellBasedEnergyRing[1], &aGroup, folderName + "_cellBasedEnergyRing1", "Cell Based Energy in Ring 2: 0.05 < R < 0.075", 50, 0, 100);
  if (sc.isSuccess()){
    folder.h_dRmax->GetXaxis()->SetTitle("Energy (GeV)");
    folder.h_dRmax->GetYaxis()->SetTitle("Number of Candidates");
  }
  */

  sc = Book2DHist ( &folder.h_centFrac_vs_LB, &aGroup, folderName + "_centFracVsLB", "Centrality Fraction vs Lumiblock", 51, 0.0, 1.02, m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0 );
  if (sc.isSuccess()){
    folder.h_centFrac_vs_LB->GetXaxis()->SetTitle("Centrality Fraction");
    folder.h_centFrac_vs_LB->GetYaxis()->SetTitle("Lumiblock");
  }

  sc = Book2DHist ( &folder.h_isolFrac_vs_LB, &aGroup, folderName + "_isolFracVsLB", "Isolation Fraction vs Lumiblock", 51, 0.0, 1.02, m_maxNLB/10+1, -5.0, (double)m_maxNLB+5.0 );
  if (sc.isSuccess()){
    folder.h_centFrac_vs_LB->GetXaxis()->SetTitle("Isolation Fraction");
    folder.h_centFrac_vs_LB->GetYaxis()->SetTitle("Lumiblock");
  }

  return sc;
}


 StatusCode tauMonTool::bookCombHistos(s_combFolder& folder,std::string folderName, Interval_t interval) {
  
  folderName = folderName + "/Combined";
  MonGroup aGroup(this, folderName, interval);
  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;

  sc = Book1DHist ( &folder.h_etOverPtLeadTrack, &aGroup, folderName + "_etOverPtLeadTrack","Et over Pt of lead track of tau candidates", 50, 0.0, 10.0);
  if (sc.isSuccess()) {
    folder.h_etOverPtLeadTrack->GetXaxis()->SetTitle("Et/Pt");
    folder.h_etOverPtLeadTrack->GetYaxis()->SetTitle("Number of Candidates");
  }

  return sc;
}

StatusCode tauMonTool::bookBDTLooseHistos(s_BDTFolder& folder, std::string folderName, Interval_t interval) {
  
  folderName = folderName + "/BDTLoose";
  MonGroup aGroup(this,folderName,interval);
  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;

  sc = Book1DHist ( &folder.h_et, &aGroup, folderName + "_et", "Et of Tau Candidate", 50, 0.0, 250.0);
  if(sc.isSuccess()){
    folder.h_et->GetXaxis()->SetTitle("Et (GeV)");
    folder.h_et->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_eta, &aGroup, folderName + "_eta", "Eta of Tau Candidate", 51, -2.55, 2.55);
  if(sc.isSuccess()){
    folder.h_eta->GetXaxis()->SetTitle("Eta");
    folder.h_eta->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_phi, &aGroup, folderName + "_phi", "Phi of Tau Candidate", 65, PHIMIN+PHIMIN/64, PHIMAX+PHIMAX/64);
  if(sc.isSuccess()){
    folder.h_phi->GetXaxis()->SetTitle("Phi");
    folder.h_phi->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_nTracks, &aGroup, folderName + "_numTracks", "Number of Tracks", 21, -0.5, 20.5);
  if(sc.isSuccess()){
    folder.h_nTracks->GetXaxis()->SetTitle("Number of Tracks");
    folder.h_nTracks->GetYaxis()->SetTitle("Number of Candidates"); 
  }

  return sc;
}

StatusCode tauMonTool::bookBDTMedHistos(s_BDTFolder& folder, std::string folderName, Interval_t interval) {
  
  folderName = folderName + "/BDTMedium";
  MonGroup aGroup(this,folderName,interval);
  folderName = this->fixName(folderName);

  StatusCode sc = StatusCode::SUCCESS;

  sc = Book1DHist ( &folder.h_et, &aGroup, folderName + "_et", "Et of Tau Candidate", 50, 0.0, 250.0);
  if(sc.isSuccess()){
    folder.h_et->GetXaxis()->SetTitle("Et (GeV)");
    folder.h_et->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_eta, &aGroup, folderName + "_eta", "Eta of Tau Candidate", 51, -2.55, 2.55);
  if(sc.isSuccess()){
    folder.h_eta->GetXaxis()->SetTitle("Eta");
    folder.h_eta->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_phi, &aGroup, folderName + "_phi", "Phi of Tau Candidate", 65, PHIMIN+PHIMIN/64, PHIMAX+PHIMAX/64);
  if(sc.isSuccess()){
    folder.h_phi->GetXaxis()->SetTitle("Phi");
    folder.h_phi->GetYaxis()->SetTitle("Number of Candidates");
  }

  sc = Book1DHist ( &folder.h_nTracks, &aGroup, folderName + "_numTracks", "Number of Tracks", 21,-0.5, 20.5);
  if(sc.isSuccess()){
    folder.h_nTracks->GetXaxis()->SetTitle("Number of Tracks");
    folder.h_nTracks->GetYaxis()->SetTitle("Number of Candidates"); 
  }

  return sc;
}

StatusCode tauMonTool::fillHistograms(s_mainFolder& mainFolder, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau, 
				      const Analysis::TauCommonDetails* p_taudetails, const Analysis::TauPID* p_tauid, 
				      const Trk::TrackSummary* trackSummary, const Trk::Perigee* perigee) {
  //MsgStream log( msgSvc(), name() );
  //log << MSG::INFO << "StatusCode tauMonTool::fillHistograms" << endreq;

  if ( fillKinHistos  (mainFolder.kinFolder, author, ftau) . isFailure()) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to fill kinematic histograms" << endreq;
  }

  if ( fillIDHistos   (mainFolder.idFolder, p_tauid, ftau) . isFailure()) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to fill identification histograms" << endreq;
  }

  if ( fillTrackHistos(mainFolder.trkFolder, p_taudetails, trackSummary, perigee) . isFailure()) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to fill traking histograms" << endreq;
  }
    
  if ( fillCaloHistos (mainFolder.caloFolder, ftau, p_taudetails) . isFailure()) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to fill calorimeter histograms" << endreq;
  }

  if ( fillCombHistos (mainFolder.combFolder, p_taudetails) . isFailure()) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to fill combined histograms" << endreq;
  }
  
  return  StatusCode::SUCCESS;
}

StatusCode tauMonTool::fillKinHistos(s_kinFolder& folder, std::vector<int> author, Analysis::TauJetContainer::const_iterator ftau) {
 
  for (unsigned int iA=0; iA < author.size(); iA++)
    folder.h_author->Fill(author[iA]);

  float eta    = (*ftau)->eta();
  float et     = (*ftau)->et() / 1000.0;
  float phi    = (*ftau)->phi();
  int   charge = (int) (*ftau)->charge();
    
  //folder.h_ntaus->    
  folder.h_eta->Fill( eta );
  folder.h_phi->Fill( phi );
  folder.h_et->Fill( et );
  folder.h_charge->Fill( charge);
  folder.h_PhiVsEta->Fill( eta, phi);
  folder.h_Eta_vs_LB->Fill( eta, m_currentLB);
  folder.h_Phi_vs_LB->Fill( phi, m_currentLB);

  return StatusCode::SUCCESS;
}


StatusCode tauMonTool::fillIDHistos(s_idFolder& folder,  const Analysis::TauPID* p_tauid, Analysis::TauJetContainer::const_iterator ftau) {
  
  folder.h_tauCutLoose         ->Fill(p_tauid->isTau(TauJetParameters::TauCutLoose));
  folder.h_tauCutMedium        ->Fill(p_tauid->isTau(TauJetParameters::TauCutMedium));
  folder.h_tauCutTight         ->Fill(p_tauid->isTau(TauJetParameters::TauCutTight));
  folder.h_tauLlhLoose         ->Fill(p_tauid->isTau(TauJetParameters::TauLlhLoose));
  folder.h_tauLlhMedium        ->Fill(p_tauid->isTau(TauJetParameters::TauLlhMedium));
  folder.h_tauLlhTight         ->Fill(p_tauid->isTau(TauJetParameters::TauLlhTight));
  folder.h_electronVetoLoose   ->Fill(p_tauid->isTau(TauJetParameters::ElectronVetoLoose));
  folder.h_electronVetoMedium  ->Fill(p_tauid->isTau(TauJetParameters::ElectronVetoMedium));
  folder.h_electronVetoTight   ->Fill(p_tauid->isTau(TauJetParameters::ElectronVetoTight));
  folder.h_muonVeto            ->Fill(p_tauid->isTau(TauJetParameters::MuonVeto));
  folder.h_eleBDTLoose         ->Fill(p_tauid->isTau(TauJetParameters::EleBDTLoose));
  folder.h_eleBDTMedium        ->Fill(p_tauid->isTau(TauJetParameters::EleBDTMedium));
  folder.h_eleBDTTight         ->Fill(p_tauid->isTau(TauJetParameters::EleBDTTight));
  folder.h_BDTJetScore         ->Fill(p_tauid->discriminant(TauJetParameters::BDTJetScore));
  folder.h_BDTEleScore         ->Fill(p_tauid->discriminant(TauJetParameters::BDTEleScore));
  folder.h_BDTJetScoreSigTrans ->Fill(p_tauid->discriminant(TauJetParameters::BDTJetScoreSigTrans));
  folder.h_BDTJetScoreSigTrans ->Fill(p_tauid->discriminant(TauJetParameters::BDTJetScoreBkgTrans));

  if(p_tauid->isTau(TauJetParameters::EleBDTLoose) == 0){
    if ( fillBDTHistos   (folder.BDTLooseFolder, ftau) . isFailure()) {
      MsgStream log( msgSvc(), name() );
      log << MSG::ERROR << "Failed to fill BDTLoose histograms" << endreq;
    }
  }
  
  if(p_tauid->isTau(TauJetParameters::EleBDTMedium) == 0){
    if ( fillBDTHistos   (folder.BDTMedFolder, ftau) . isFailure()) {
      MsgStream log( msgSvc(), name() );
      log << MSG::ERROR << "Failed to fill BDTMed histograms" << endreq;
    }
  }

  return StatusCode::SUCCESS;
 
}

StatusCode tauMonTool::fillTrackHistos(s_trkFolder& folder, const Analysis::TauCommonDetails* p_taudetails, const Trk::TrackSummary* trackSummary, 
				       const Trk::Perigee* perigee) {

  if( p_taudetails ) {
    folder.h_massTrkSys            ->Fill( p_taudetails->massTrkSys() / 1000.0); //in GeV
    folder.h_trkWidth2             ->Fill( p_taudetails->trkWidth2() );
    folder.h_trFlightPathSig       ->Fill( p_taudetails->trFlightPathSig() );
    folder.h_ipSigLeadTrk          ->Fill( p_taudetails->ipSigLeadTrk() );
    folder.h_ipZ0SinThetaSigLeadTrk->Fill( p_taudetails->ipZ0SinThetaSigLeadTrk() );
    folder.h_leadTrkPt             ->Fill( p_taudetails->leadTrkPt()/ 1000.0 );
  }
  
  if ( trackSummary ) {    
    folder.h_numberOfTRTHighThresholdHits    ->Fill(trackSummary->get( Trk::numberOfTRTHighThresholdHits) );
    folder.h_numberOfTRTHits                 ->Fill(trackSummary->get( Trk::numberOfTRTHits ));
    folder.h_numberOfSCTHits                 ->Fill(trackSummary->get( Trk::numberOfSCTHits ));	
    folder.h_numberOfPixelHits               ->Fill(trackSummary->get( Trk::numberOfPixelHits ));
    folder.h_numberOfTRTHighThresholdOutliers->Fill(trackSummary->get( Trk::numberOfTRTHighThresholdOutliers ));
    folder.h_numberOfTRTOutliers             ->Fill(trackSummary->get( Trk::numberOfTRTOutliers ));
  }
  
  if (perigee) {
    folder.h_d0                              ->Fill(perigee->parameters()[Trk::d0]);
    folder.h_z0                              ->Fill(perigee->parameters()[Trk::z0]);
    folder.h_phi                             ->Fill(perigee->parameters()[Trk::phi]);
    folder.h_eta                             ->Fill(perigee->eta());
    folder.h_pT                              ->Fill(perigee->pT()/ 1000.0);
    folder.h_z0_vs_LB                        ->Fill(perigee->parameters()[Trk::z0], m_currentLB);
  }

  return StatusCode::SUCCESS;
}

StatusCode tauMonTool::fillCaloHistos(s_caloFolder& folder, Analysis::TauJetContainer::const_iterator ftau, const Analysis::TauCommonDetails* p_taudetails) {

  folder.h_eta                   ->Fill(( *ftau)->eta());
  folder.h_phi                   ->Fill((*ftau)->phi());
  folder.h_etEMCalib             ->Fill( p_taudetails->seedCalo_etEMCalib() / 1000.0);
  folder.h_etHadCalib            ->Fill( p_taudetails->seedCalo_etHadCalib() / 1000.0);
  folder.h_etEMAtEMScale         ->Fill( p_taudetails->seedTrk_etEMAtEMScale() / 1000.0);
  folder.h_etHadAtEMScale        ->Fill( p_taudetails->seedTrk_etHadAtEMScale() / 1000.0);
  folder.h_EMRadius              ->Fill( p_taudetails->seedCalo_EMRadius());
  folder.h_centFrac              ->Fill( p_taudetails->seedCalo_centFrac());
  folder.h_hadRadius             ->Fill( p_taudetails->seedCalo_hadRadius());
  folder.h_isolFrac              ->Fill( p_taudetails->seedCalo_isolFrac());
  folder.h_stripWidth2           ->Fill( p_taudetails->seedCalo_stripWidth2());
  folder.h_nStrip                ->Fill( p_taudetails->seedCalo_nStrip());
  folder.h_trkAvgDist            ->Fill( p_taudetails->seedCalo_trkAvgDist());
  folder.h_dRmax                 ->Fill( p_taudetails->seedCalo_dRmax());
  /*
  folder.h_cellBasedEnergyRing[0]->Fill( p_taudetails->cellBasedEnergyRing1());
  folder.h_cellBasedEnergyRing[1]->Fill( p_taudetails->cellBasedEnergyRing2());
  */
  folder.h_centFrac_vs_LB        ->Fill( p_taudetails->seedCalo_centFrac(), m_currentLB);
  folder.h_isolFrac_vs_LB        ->Fill( p_taudetails->seedCalo_isolFrac(), m_currentLB);

  return  StatusCode::SUCCESS;
}


StatusCode tauMonTool::fillCombHistos(s_combFolder& folder, const Analysis::TauCommonDetails* p_taudetails) {
    
  folder.h_etOverPtLeadTrack->Fill(p_taudetails->etOverPtLeadTrk());

  return StatusCode::SUCCESS;
}

StatusCode tauMonTool::fillBDTHistos(s_BDTFolder& folder, Analysis::TauJetContainer::const_iterator ftau) {
  float eta     = (*ftau)->eta();
  float et      = (*ftau)->et() / 1000.0;
  float phi     = (*ftau)->phi();
  float nTracks = (*ftau)->numTrack();

  folder.h_et     ->Fill(et);
  folder.h_eta    ->Fill(eta);
  folder.h_phi    ->Fill(phi);
  folder.h_nTracks->Fill(nTracks);

  return StatusCode::SUCCESS;
}


StatusCode tauMonTool::fillPhysicsHistograms(Analysis::TauJetContainer::const_iterator ftau, const Trk::Perigee* perigee){

  float eta = (*ftau)->eta();
  
  m_eta_Tau_Z->Fill(eta);

  if(perigee){
    m_pTVsEta_Tau_Z->Fill(perigee->pT()/1000.0, eta);
    m_pTVsEta_Tau_W->Fill(perigee->pT()/1000.0, eta);
  }

  return StatusCode::SUCCESS;

}


StatusCode tauMonTool::bookPhysicsHistograms(){
  MsgStream log(msgSvc(), name() );
  log << MSG::INFO << "StatusCode tauMonTool::bookPhysicsHistograms" << endreq;

  StatusCode sc = StatusCode::SUCCESS;

  //**********
  // Z details
  //**********
  MonGroup tau_Z (this, "Tau/Physics/Z", run);

  sc = Book1DHist (&m_eta_Tau_Z, &tau_Z, "tau_eta", "Eta of Tau Candidates", 51, -2.55, 2.55);
  m_eta_Tau_Z->GetXaxis()->SetTitle("Eta");
  sc = Book2DHist (&m_pTVsEta_Tau_Z, &tau_Z, "tau_pTVsEta", "Pt vs. Eta of Tau Candidates", 51, -2.55, 2.55, 100, 0.0, 150.0);
  m_pTVsEta_Tau_Z->GetXaxis()->SetTitle("Eta");
  m_pTVsEta_Tau_Z->GetYaxis()->SetTitle("Pt (GeV)");
  sc = Book2DHist (&m_pTVsEta_Lepton_Z, &tau_Z, "lepton_pTVsEta", "Pt vs. Eta of Lepton Candidates", 51, -2.55, 2.55, 100, 0.0, 150.0);
  m_pTVsEta_Lepton_Z->GetXaxis()->SetTitle("Eta");
  m_pTVsEta_Lepton_Z->GetYaxis()->SetTitle("Pt (GeV)");

  //**********
  // W details
  //**********
  MonGroup tau_W (this, "Tau/Physics/W", run);

  sc = Book2DHist (&m_pTVsEta_Tau_W, &tau_W, "tau_pTVsEta", "Pt vs. Eta of Tau Candidates", 51, -2.55, 2.55, 100, 0.0, 150.0);
  m_pTVsEta_Tau_W->GetXaxis()->SetTitle("Eta");
  m_pTVsEta_Tau_W->GetYaxis()->SetTitle("Pt (GeV)");

  return StatusCode::SUCCESS;
}


/*
StatusCode tauMonTool::fillOldHistograms(Analysis::TauJetContainer::const_iterator ftau, const Analysis::TauCommonDetails* p_taudetails, const Analysis::TauPID* p_tauid) {
  

  //--------------------
  // To keep track of type, tauRec-only, 1p3p-only, or merged
  // merged       -> auth_Case = 0
  // tauRec-only  -> auth_Case = 1
  // 1p3p-only    -> auth_Case = 2 
  //--------------------
  int auth_Case = -1;
  bool isCalo = (*ftau)->hasAuthor( TauJetParameters::tauRec);
  bool isTrk  = (*ftau)->hasAuthor( TauJetParameters::tau1P3P);

  if (isCalo && isTrk)  auth_Case = 0;
  if (isCalo && !isTrk) auth_Case = 1;
  if (!isCalo && isTrk) auth_Case = 2;


  //eta, phi
  float eta    = (*ftau)->eta();
  float phi    = (*ftau)->phi();

   
  // merged category (both track and calo seeds available)
  
  if (auth_Case == 0) // tauRec & tau1p3p
    {
      
      //Access the ID flags:
            
      m_tauCutLoose_MERG->Fill(p_tauid->isTau(TauJetParameters::TauCutLoose));
      m_tauCutMedium_MERG->Fill(p_tauid->isTau(TauJetParameters::TauCutMedium));
      m_tauCutTight_MERG->Fill(p_tauid->isTau(TauJetParameters::TauCutTight));
      m_electronVeto_MERG->Fill(p_tauid->isTau(TauJetParameters::ElectronVeto));
      m_muonVeto_MERG->Fill(p_tauid->isTau(TauJetParameters::MuonVeto));
      m_tauLlhTight_MERG->Fill(p_tauid->isTau(TauJetParameters::TauLlhTight));
      m_tauLlhMedium_MERG->Fill(p_tauid->isTau(TauJetParameters::TauLlhMedium));
      m_tauLlhLoose_MERG->Fill(p_tauid->isTau(TauJetParameters::TauLlhLoose));
      
      // Access tauCommon details
      //const Analysis::TauCommonDetails*  p_taudetails = (*ftau)->details<const Analysis::TauCommonDetails>();	    
      
      //const Analysis::TauCommonDetails*  p_taudetails_1P3P = (*ftau)->details<const Analysis::TauCommonDetails>();
      
      
      double phiCalo_MERG = p_taudetails->seedCalo_phi();
      m_phiCalo_MERG->Fill(phiCalo_MERG);
      
      double etaCalo_MERG = p_taudetails->seedCalo_eta();
      m_etaCalo_MERG->Fill(etaCalo_MERG);
      
      double emRadius_MERG = p_taudetails->seedCalo_EMRadius(); 
      m_emRadius_MERG->Fill(emRadius_MERG);
	    
      double isolationFraction_MERG  = p_taudetails->seedCalo_isolFrac();
      m_isolationFraction_MERG->Fill(isolationFraction_MERG);

      double centralityFraction_MERG = p_taudetails->seedCalo_centFrac();
      m_centralityFraction_MERG->Fill(centralityFraction_MERG);

      double stripWidth2_MERG = p_taudetails->seedCalo_stripWidth2();
      m_stripWidth2_MERG->Fill(stripWidth2_MERG);
	    
      double numStripCells_MERG = p_taudetails->seedCalo_nStrip();
      m_numStripCells_MERG->Fill(numStripCells_MERG);

      double etEMCalib_MERG = p_taudetails->seedCalo_etEMCalib();
      m_etEMCalib_MERG->Fill(etEMCalib_MERG / 1000.0); // in GeV

      double etHadCalib_MERG = p_taudetails->seedCalo_etHadCalib();
      m_etHadCalib_MERG->Fill(etHadCalib_MERG / 1000.0 ); // in GeV

      double etOverPtLeadTrack_MERG = p_taudetails->etOverPtLeadTrk();
      m_etOverPtLeadTrack_MERG->Fill(etOverPtLeadTrack_MERG);

      double trFlightPathSig_MERG = p_taudetails->trFlightPathSig();
      m_trFlightPathSig_MERG->Fill(trFlightPathSig_MERG);

      double leadingTrackPT_MERG =  p_taudetails->leadTrkPt();
      m_leadingTrackPT_MERG->Fill(leadingTrackPT_MERG / 1000.0 ); // in GeV

      double ipSigLeadTrack_MERG = p_taudetails->ipSigLeadTrk();
      m_ipSigLeadTrack_MERG->Fill(ipSigLeadTrack_MERG);

      //2-D histograms:
      m_EtOverPtLTVsEta_MERG->Fill(etaCalo_MERG,etOverPtLeadTrack_MERG);
      m_EtOverPtLTVsPhi_MERG->Fill(phiCalo_MERG,etOverPtLeadTrack_MERG);
      m_etEmVsEta_MERG->Fill(etaCalo_MERG,etEMCalib_MERG / 1000.0);
      m_etHadVsEta_MERG->Fill(etaCalo_MERG,etHadCalib_MERG / 1000.0);
      m_etEmVsPhi_MERG->Fill(phiCalo_MERG,etEMCalib_MERG / 1000.0);
      m_etHadVsPhi_MERG->Fill(phiCalo_MERG,etHadCalib_MERG / 1000.0);

      double etChrgHAD_MERG = p_taudetails->seedTrk_etChrgHad();
      m_etChrgHAD_MERG->Fill(etChrgHAD_MERG / 1000.0 ); // in GeV

      double etIsolEM_MERG = p_taudetails->seedTrk_etIsolEM();
      m_etIsolEM_MERG->Fill(etIsolEM_MERG / 1000.0 );  // in Gev

      double etIsolHAD_MERG = p_taudetails->seedTrk_etIsolHad();
      m_etIsolHAD_MERG->Fill(etIsolHAD_MERG / 1000.0 );  // in GeV

      int nAssocTracksCore_MERG = p_taudetails->seedTrk_nOtherCoreTrk();
      m_nAssocTracksCore_MERG->Fill(nAssocTracksCore_MERG);

      int nAssocTracksIsol_MERG = p_taudetails->seedTrk_nIsolTrk();
      m_nAssocTracksIsol_MERG->Fill(nAssocTracksIsol_MERG);

      double massTrk3P_MERG = p_taudetails->massTrkSys();
      m_massTrk3P_MERG->Fill(massTrk3P_MERG / 1000.0 );  // in Gev

      double rWidth2Trk3P_MERG = p_taudetails->trkWidth2();
      m_rWidth2Trk3P_MERG->Fill(rWidth2Trk3P_MERG);

      double signD0Trk3P_MERG = p_taudetails->ipSigLeadTrk();
      m_signD0Trk3P_MERG->Fill(signD0Trk3P_MERG);

      double etHadAtEMScale_MERG = p_taudetails->seedTrk_etHadAtEMScale();
      m_etHadAtEMScale_MERG->Fill(etHadAtEMScale_MERG / 1000.0 );  // in GeV

      double etEMAtEMScale_MERG = p_taudetails->seedTrk_etEMAtEMScale();
      m_etEMAtEMScale_MERG->Fill(etEMAtEMScale_MERG / 1000.0 );  // in GeV

      double etEMCL_MERG = p_taudetails->seedTrk_etEMCL();
      m_etEMCL_MERG->Fill(etEMCL_MERG / 1000.0 ); // in GeV

      double etChrgEM_MERG = p_taudetails->seedTrk_etChrgEM();
      m_etChrgEM_MERG->Fill(etChrgEM_MERG / 1000.0 );  // in GeV

      double etNeuEM_MERG = p_taudetails->seedTrk_etNeuEM();
      m_etNeuEM_MERG->Fill(etNeuEM_MERG / 1000.0 );  // in GeV

      double etResNeuEM_MERG = p_taudetails->seedTrk_etResNeuEM();
      m_etResNeuEM_MERG->Fill(etResNeuEM_MERG);

      double numPi0_MERG = p_taudetails->nPi0();
      m_numPi0_MERG->Fill(numPi0_MERG);

      double etChrgHADoverPttot_MERG = p_taudetails->seedTrk_etChrgHadOverSumTrkPt();
      m_etChrgHADoverPttot_MERG->Fill(etChrgHADoverPttot_MERG);

      double z0SinThetaSig_MERG = p_taudetails->ipZ0SinThetaSigLeadTrk();
      m_z0SinThetaSig_MERG->Fill(z0SinThetaSig_MERG);

      double etIsolFrac_MERG = p_taudetails->seedTrk_isolFracWide();
      m_etIsolFrac_MERG->Fill(etIsolFrac_MERG);

      // Details regarding tracks (1P3P-only)

      for(unsigned int n = 1; n <= 3; n++){
	if((*ftau)->numTrack() >= n)  
	  {
	    const Trk::TrackSummary *summary = (*ftau)->track(n-1)->trackSummary();
		  
	    if ( !summary)
	      {
		return StatusCode::FAILURE;
	      }
		  
	    else
	      {
		int trtHighHits_MERG = -1;
		trtHighHits_MERG = summary->get( Trk::numberOfTRTHighThresholdHits );
		m_HighTRTHits_MERG->Fill(trtHighHits_MERG);	
		      
		int trtLowHits_MERG = -1;
		trtLowHits_MERG = summary->get( Trk::numberOfTRTHits );
		m_LowTRTHits_MERG->Fill(trtLowHits_MERG);	
		      
		int SCTHits_MERG = -1;
		SCTHits_MERG = summary->get( Trk::numberOfSCTHits );
		m_SCTHits_MERG->Fill(SCTHits_MERG);	
		      
		int PixelHits_MERG = -1;
		PixelHits_MERG = summary->get( Trk::numberOfPixelHits );
		m_PixelHits_MERG->Fill(PixelHits_MERG);	    
	      }
	  }  
      } // have at least 1,2 or 3 tracks	    

      //2-D Plots
	    
      m_numPi0VsEta_MERG->Fill(eta,numPi0_MERG);
      m_numPi0VsPhi_MERG->Fill(phi,numPi0_MERG);
      m_nTracksCoreVsEta_MERG->Fill(eta,nAssocTracksCore_MERG);
      m_nTracksCoreVsPhi_MERG->Fill(phi,nAssocTracksCore_MERG);
      m_etIsolFracVsEta_MERG->Fill(eta,etIsolFrac_MERG ); 	    
      m_etIsolFracVsPhi_MERG->Fill(phi,etIsolFrac_MERG );
      
      
    }

  return StatusCode::SUCCESS;
}




StatusCode tauMonTool::bookOldHistograms() {
  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "StatusCode tauMonTool::bookOldHistograms" << endreq;

  StatusCode sc = StatusCode::SUCCESS;

      //-------------------
      //- Merged Candidates
      //-------------------
      MonGroup tau_MERG (this, "Tau/Expert/Merged", expert, run);      
      sc = Book1DHist ( &m_etaCalo_MERG, &tau_MERG, "MERG_etaCalo","Calorimeter eta of tau candidates", 51, -2.55, 2.55);
      m_etaCalo_MERG->GetXaxis()->SetTitle("Eta");
      m_etaCalo_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_phiCalo_MERG, &tau_MERG, "MERG_phiCalo","Calorimeter phi of tau candidates", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64.);
      m_phiCalo_MERG->GetXaxis()->SetTitle("Phi");
      m_phiCalo_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_emRadius_MERG, &tau_MERG, "MERG_emRadius","EM Radius of tau candidates", 50, 0.0, 1.0);
      m_emRadius_MERG->GetXaxis()->SetTitle("EM radius");
      m_emRadius_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_isolationFraction_MERG, &tau_MERG, "MERG_isolationFraction","Isolation fraction of tau candidates", 50, 0.0, 1.0);
      m_isolationFraction_MERG->GetXaxis()->SetTitle("Isolation Fraction");
      m_isolationFraction_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_centralityFraction_MERG, &tau_MERG, "MERG_centralityFraction","Centrality fraction of tau candidates", 50, 0.0, 1.0);
      m_centralityFraction_MERG->GetXaxis()->SetTitle("Centrality Fraction");
      m_centralityFraction_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_stripWidth2_MERG, &tau_MERG, "MERG_stripWidth2","Strip Width of tau candidates", 50, -0.1, 0.1);
      m_stripWidth2_MERG->GetXaxis()->SetTitle("Strip Width");
      m_stripWidth2_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_numStripCells_MERG, &tau_MERG, "MERG_numStripCells","Number of strip cells", 50, -0.5, 49.5);
      m_numStripCells_MERG->GetXaxis()->SetTitle("Number of Strip Cells");
      m_numStripCells_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etEMCalib_MERG, &tau_MERG, "MERG_etEMCalib","Calibrated EM Et of tau candidates", 50, 0.0, 150.0);
      m_etEMCalib_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etEMCalib_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etHadCalib_MERG, &tau_MERG, "MERG_etHadCalib","Calibrated Had Et of tau candidates", 50, 0.0, 150.0);
      m_etHadCalib_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etHadCalib_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etOverPtLeadTrack_MERG, &tau_MERG, "MERG_etOverPtLeadTrack","Et over Pt of lead track of tau candidates", 50, 0.0, 10.0);
      m_etOverPtLeadTrack_MERG->GetXaxis()->SetTitle("Et/Pt");
      m_etOverPtLeadTrack_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_trFlightPathSig_MERG, &tau_MERG, "MERG_trFlightPathSig","Transverse flight path of tau candidates", 50, -5.0, 5.0);
      m_trFlightPathSig_MERG->GetXaxis()->SetTitle("Flight Path Significance");
      m_trFlightPathSig_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_leadingTrackPT_MERG, &tau_MERG, "MERG_leadingTrackPT","Pt of leading track of tau candidates", 50, 0.0, 150.0);
      m_leadingTrackPT_MERG->GetXaxis()->SetTitle("Pt (GeV)");
      m_leadingTrackPT_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_ipSigLeadTrack_MERG, &tau_MERG, "MERG_ipSigLeadTrack","IP Significance of leading track of tau candidates", 50, -2.0, 5.0);
      m_ipSigLeadTrack_MERG->GetXaxis()->SetTitle("IP Significance");
      m_ipSigLeadTrack_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_tauCutLoose_MERG, &tau_MERG, "MERG_tauCutLoose","ID Flag: TauCutLoose", 2, -0.5, 1.5);
      m_tauCutLoose_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauCutLoose_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_tauCutMedium_MERG, &tau_MERG, "MERG_tauCutMedium","ID Flag: TauCutMedium", 2, -0.5, 1.5);
      m_tauCutMedium_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauCutMedium_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_tauCutTight_MERG, &tau_MERG, "MERG_tauCutTight","ID Flag: TauCutTight", 2, -0.5, 1.5);
      m_tauCutTight_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauCutTight_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_electronVeto_MERG, &tau_MERG, "MERG_electronVeto","ID Flag: Electron Veto", 2, -0.5, 1.5);
      m_electronVeto_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_electronVeto_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_muonVeto_MERG, &tau_MERG, "MERG_muonVeto","ID Flag: Muon Veto", 2, -0.5, 1.5);
      m_muonVeto_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_muonVeto_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_tauLlhTight_MERG, &tau_MERG, "MERG_tauLlhTight","ID Flag: Tau Likelihood Tight", 2, -0.5, 1.5);
      m_tauLlhTight_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauLlhTight_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_tauLlhMedium_MERG, &tau_MERG, "MERG_tauLlhMedium","ID Flag: Tau Likelihood Medium", 2, -0.5, 1.5);
      m_tauLlhMedium_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauLlhMedium_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_tauLlhLoose_MERG, &tau_MERG, "MERG_tauLlhLoose","ID Flag: Tau Likelihood Loose", 2, -0.5, 1.5);
      m_tauLlhLoose_MERG->GetXaxis()->SetBinLabel(1, "False");
      m_tauLlhLoose_MERG->GetXaxis()->SetBinLabel(2, "True");
      sc = Book1DHist ( &m_etChrgHAD_MERG, &tau_MERG, "MERG_etChrgHAD","Charged Had Et of tau candidates", 50, 0.0, 150.0);
      m_etChrgHAD_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etChrgHAD_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etIsolEM_MERG, &tau_MERG, "MERG_etIsolEM","Et in EM Calo in Isolation Region", 50, 0.0, 150.0);
      m_etIsolEM_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etIsolEM_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etIsolHAD_MERG, &tau_MERG, "MERG_etIsolHAD","Et in HAD Calo in Isolation Region", 50, 0.0, 150.0);
      m_etIsolHAD_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etIsolHAD_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_nAssocTracksCore_MERG, &tau_MERG, "MERG_nAssocTracksCore","Number of Associated Tracks before Quality Cuts", 16, -0.5, 15.5);
      m_nAssocTracksCore_MERG->GetXaxis()->SetTitle("Number of Tracks");
      m_nAssocTracksCore_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_nAssocTracksIsol_MERG, &tau_MERG, "MERG_nAssocTracksIsol","Number of Tracks in Isolation Region", 16, -0.5, 15.5);
      m_nAssocTracksIsol_MERG->GetXaxis()->SetTitle("Number of Tracks");
      m_nAssocTracksIsol_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_massTrk3P_MERG, &tau_MERG, "MERG_massTrk3P","Invariant mass of the track system", 30, -1.0, 5.0);
      m_massTrk3P_MERG->GetXaxis()->SetTitle("Mass (GeV)");
      m_massTrk3P_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_rWidth2Trk3P_MERG, &tau_MERG, "MERG_rWidth2Trk3P","Width of track momenta", 50, 0.0, 1.0);
      m_rWidth2Trk3P_MERG->GetXaxis()->SetTitle("Momenta (GeV)");
      m_rWidth2Trk3P_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_signD0Trk3P_MERG, &tau_MERG, "MERG_signD0Trk3P","Signed Transverse Impact Parameter", 50, -5.0, 15.0);
      m_signD0Trk3P_MERG->GetXaxis()->SetTitle("Impact Paraeter (mm)");
      m_signD0Trk3P_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etHadAtEMScale_MERG, &tau_MERG, "MERG_etHadAtEMScale","Hadronic energy at the EM Scale", 50, 0.0, 150.0);
      m_etHadAtEMScale_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etHadAtEMScale_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etEMAtEMScale_MERG, &tau_MERG, "MERG_etEMAtEMScale","EM energy at the EM scale", 50, 0.0, 150.0);
      m_etEMAtEMScale_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etEMAtEMScale_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etEMCL_MERG, &tau_MERG, "MERG_etEMCL","Et of cells classified as pure electromagnetic seeded by topo cluster", 50, 0.0, 150.0);
      m_etEMCL_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etEMCL_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etChrgEM_MERG, &tau_MERG, "MERG_etChrgEM","Charged Et of EM cells near track", 50, 0.0, 150.0);
      m_etChrgEM_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etChrgEM_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etNeuEM_MERG, &tau_MERG, "MERG_etNeuEM","Neutral Et of EM cells after subtraction", 50, 0.0, 150.0);
      m_etNeuEM_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etNeuEM_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etResNeuEM_MERG, &tau_MERG, "MERG_etResNeuEM","Correction term for eflow calculation", 50, 0.0, 150.0);
      m_etResNeuEM_MERG->GetXaxis()->SetTitle("Et (GeV)");
      m_etResNeuEM_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_numPi0_MERG, &tau_MERG, "MERG_numPi0","Number of Reconstructed Pizero Clusters", 11, -0.5, 10.5);
      m_numPi0_MERG->GetXaxis()->SetTitle("Number of Clusters");
      m_numPi0_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etChrgHADoverPttot_MERG, &tau_MERG, "MERG_etChrgHADoverPttot","Charged HAD Et over sum of track Pt", 50, 0.0, 150.0);
      m_etChrgHADoverPttot_MERG->GetXaxis()->SetTitle("Et/Pt");
      m_etChrgHADoverPttot_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_z0SinThetaSig_MERG, &tau_MERG, "MERG_z0SinThetaSig","Significance of z0SinTheta", 50, -10.0, 10.0);
      m_z0SinThetaSig_MERG->GetXaxis()->SetTitle("Significance");
      m_z0SinThetaSig_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_etIsolFrac_MERG, &tau_MERG, "MERG_etIsolFrac","Et Isolation Fraction", 50, 0.0, 4.0);
      m_etIsolFrac_MERG->GetXaxis()->SetTitle("Isolation Fraction");
      m_etIsolFrac_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_HighTRTHits_MERG, &tau_MERG, "MERG_HighTRTHits","Number of High TRT Hits", 15, -0.5, 14.5);
      m_HighTRTHits_MERG->GetXaxis()->SetTitle("Number of Hits");
      m_HighTRTHits_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_LowTRTHits_MERG, &tau_MERG, "MERG_LowTRTHits","Number of Low TRT Hits", 101, -0.5, 100.5);
      m_LowTRTHits_MERG->GetXaxis()->SetTitle("Number of Hits");
      m_LowTRTHits_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_SCTHits_MERG, &tau_MERG, "MERG_SCTHits","Number of SCT Hits", 51, -0.5, 50.5);
      m_SCTHits_MERG->GetXaxis()->SetTitle("Number of Hits");
      m_SCTHits_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book1DHist ( &m_PixelHits_MERG, &tau_MERG, "MERG_PixelHits","Number of Pixel Hits", 26, -0.5, 25.5);
      m_PixelHits_MERG->GetXaxis()->SetTitle("Number of Hits");
      m_PixelHits_MERG->GetYaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_EtOverPtLTVsEta_MERG, &tau_MERG, "MERG_EtOverPtLTVsEta","Et over Pt of leading track vs. eta", 51, -2.55, 2.55, 100, 0.0, 10.0);
      m_EtOverPtLTVsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_EtOverPtLTVsEta_MERG->GetYaxis()->SetTitle("Et/Pt");
      m_EtOverPtLTVsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_EtOverPtLTVsPhi_MERG, &tau_MERG, "MERG_EtOverPtLTVsPhi","Et over Pt of leading track vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 100, 0.0, 10.0);
      m_EtOverPtLTVsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_EtOverPtLTVsPhi_MERG->GetYaxis()->SetTitle("Et/Pt");
      m_EtOverPtLTVsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etEmVsEta_MERG, &tau_MERG, "MERG_etEmVsEta","Et EM calib vs. eta", 51, -2.55, 2.55, 100, 0.0, 150.0);
      m_etEmVsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_etEmVsEta_MERG->GetYaxis()->SetTitle("Et (GeV)");
      m_etEmVsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etEmVsPhi_MERG, &tau_MERG, "MERG_etEmVsPhi","Et EM calib vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 100, 0.0, 150.0);
      m_etEmVsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_etEmVsPhi_MERG->GetYaxis()->SetTitle("Et (GeV)");
      m_etEmVsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etHadVsEta_MERG, &tau_MERG, "MERG_etHadVsEta","Et Had calib vs. eta", 51, -2.55, 2.55, 100, 0.0, 150.0);
      m_etHadVsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_etHadVsEta_MERG->GetYaxis()->SetTitle("Et (GeV)");
      m_etHadVsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etHadVsPhi_MERG, &tau_MERG, "MERG_etHadVsPhi","Et Had calib vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 100, 0.0, 150.0);
      m_etHadVsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_etHadVsPhi_MERG->GetYaxis()->SetTitle("Et (GeV)");
      m_etHadVsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_numPi0VsEta_MERG, &tau_MERG, "MERG_numPi0VsEta","Number of Reconstructed PiZero Clusters vs. eta", 51, -2.55, 2.55, 11, -0.5, 10.5);
      m_numPi0VsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_numPi0VsEta_MERG->GetYaxis()->SetTitle("Number of Pi0");
      m_numPi0VsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_numPi0VsPhi_MERG, &tau_MERG, "MERG_numPi0VsPhi","Number of Reconstructed PiZero Clusters vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 11, -0.5, 10.5);
      m_numPi0VsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_numPi0VsPhi_MERG->GetYaxis()->SetTitle("Number of Pi0");
      m_numPi0VsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_nTracksCoreVsEta_MERG, &tau_MERG, "MERG_nTracksCoreVsEta","Number of Associated Tracks bfore Quality Cuts vs. eta", 51, -2.55, 2.55, 16, -0.5, 15.5);
      m_nTracksCoreVsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_nTracksCoreVsEta_MERG->GetYaxis()->SetTitle("Number of Tracks");
      m_nTracksCoreVsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_nTracksCoreVsPhi_MERG, &tau_MERG, "MERG_nTracksCoreVsPhi","Number of Associated Tracks before Quality cuts vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 16, -0.5, 15.5 );
      m_nTracksCoreVsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_nTracksCoreVsPhi_MERG->GetYaxis()->SetTitle("Number of Tracks");
      m_nTracksCoreVsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etIsolFracVsEta_MERG, &tau_MERG, "MERG_etIsolFracVsEta","Et Isolation Fraction vs. Eta", 51, -2.55, 2.55, 50, 0.0, 4.0);
      m_etIsolFracVsEta_MERG->GetXaxis()->SetTitle("Eta");
      m_etIsolFracVsEta_MERG->GetYaxis()->SetTitle("Isolation Fraction");
      m_etIsolFracVsEta_MERG->GetZaxis()->SetTitle("Number of Candidates");
      sc = Book2DHist ( &m_etIsolFracVsPhi_MERG, &tau_MERG, "MERG_etIsolFracVsPhi","Et Isolation Fraction vs. phi", 65, PHIMIN+PHIMIN/64., PHIMAX+PHIMAX/64., 50, 0.0, 4.0);
      m_etIsolFracVsPhi_MERG->GetXaxis()->SetTitle("Phi");
      m_etIsolFracVsPhi_MERG->GetYaxis()->SetTitle("Isolation Fraction");
      m_etIsolFracVsPhi_MERG->GetZaxis()->SetTitle("Number of Candidates");

      return StatusCode::SUCCESS;
}*/
