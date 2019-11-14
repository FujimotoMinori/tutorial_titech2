// MyPackage includes
#include "MyPackageAlg.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "PATInterfaces/CorrectionCode.h" // to check the return correction code status of tools 
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/ElectronxAODHelpers.h"
#include "PhotonVertexSelection/PhotonVertexHelpers.h"



MyPackageAlg::MyPackageAlg( const std::string& name, ISvcLocator* pSvcLocator ) : AthAnalysisAlgorithm( name, pSvcLocator ),
	//declareProperty( "Property", m_nProperty = 0, "My Example Integer Property" ); //example property declaration
	m_grl ("GoodRunsListSelectionTool/grl", this),
	m_muonSelection ("CP::MuonSelectionTool", this),
	m_muonCalibrationAndSmearingTool ("CP::MuonCalibrationAndSmearingTool/MuonCorrectionTool",this),
	m_isoTool("CP::IsolationSelectionTool/IsolationSelectionTool"), 
	m_btagSelectionTool ("BTaggingSelectionTool", this),
	m_EgammaCalibrationAndSmearingTool("CP::EgammaCalibrationAndSmearingTool/EgammaCalibrationAndSmearingTool", this),
	m_fudgeMCTool("ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool",this), 
	m_photonLooseIsEMSelector("AsgPhotonIsEMSelector/PhotonSelIsEM_Loose",this),
	m_photonTightIsEMSelector("AsgPhotonIsEMSelector/PhotonSelIsEM_Tight",this),
	m_deadHVTool("AsgDeadHVCellRemovalTool/deadHVTool",this),
	m_vertexTool("CP::PhotonVertexSelectionTool/vxtTool"),
	m_pointTool("CP::PhotonPointingTool/pointingTool",this),
	m_indet_tracksel("InDet::InDetTrackSelectionTool/trackSel"),
	m_track_iso_tool( "xAOD::TrackIsolationTool/TrackIsolationTool"),
	m_isoCorrTool("CP::IsolationCorrectionTool/IsoCorrTool",this),
	m_photonFCLooseIsoTool("CP::IsolationSelectionTool/PhotonIsolationSelectionTool", this){

	}


MyPackageAlg::~MyPackageAlg() {}


StatusCode MyPackageAlg::initialize() {
	ATH_MSG_INFO ("Initializing " << name() << "...");
	//
	//This is called once, before the start of the event loop
	//Retrieves of tools you have configured in the joboptions go here
	//

	//HERE IS AN EXAMPLE
	//We will create a histogram and a ttree and register them to the histsvc
	//Remember to configure the histsvc stream in the joboptions
	//
	//m_myHist = new TH1D("myHist","myHist",100,0,100);
	//CHECK( histSvc()->regHist("/MYSTREAM/myHist", m_myHist) ); //registers histogram to output stream
	//m_myTree = new TTree("myTree","myTree");
	//CHECK( histSvc()->regTree("/MYSTREAM/SubDirectory/myTree", m_myTree) ); //registers tree to output stream inside a sub-directory

	m_tdt.setTypeAndName("Trig::TrigDecisionTool/TrigDecisionTool");
	CHECK( m_tdt.initialize() );

	std::vector<std::string> vecgrl;
	vecgrl.push_back("GoodRunsLists/data15_13TeV/20170619/physics_25ns_21.0.19.xml");
	vecgrl.push_back("GoodRunsLists/data16_13TeV/20180129/physics_25ns_21.0.19.xml");
	vecgrl.push_back("GoodRunsLists/data17_13TeV/20180619/physics_25ns_Triggerno17e33prim.xml");
	vecgrl.push_back("GoodRunsLists/data18_13TeV/20190318/physics_25ns_Triggerno17e33prim.xml");
	CHECK (m_grl.setProperty("GoodRunsListVec", vecgrl));
	CHECK (m_grl.setProperty("PassThrough", false));
	CHECK (m_grl.initialize());

	CHECK (m_muonSelection.setProperty("MaxEta", 2.5));
	CHECK (m_muonSelection.setProperty("MuQuality", (int) xAOD::Muon::Quality::Medium));
	CHECK (m_muonSelection.initialize());
	CHECK (m_muonCalibrationAndSmearingTool.initialize());

	CHECK( m_isoTool.setProperty("MuonWP", "FixedCutTightTrackOnly") );
	CHECK( m_isoTool.retrieve() );

    //histograms
	m_histAverageIntPerXing = new TH1D("AverageIntPerXing","AverageIntPerXing",100,0,100);
	CHECK( histSvc()->regHist("/MYSTREAM/AverageIntPerXing", m_histAverageIntPerXing) );
	m_histdiMuon = new TH1D("diMuon","diMuon",150,0,150000);
	CHECK( histSvc()->regHist("/MYSTREAM/diMuon", m_histdiMuon) );
	m_histjet = new TH1D("jet","jet",150,0,150000);
	CHECK( histSvc()->regHist("/MYSTREAM/jet", m_histjet) );
	m_histjetbtagged = new TH1D("jetbtagged","jetbtagged",150,0,150000);
	CHECK( histSvc()->regHist("/MYSTREAM/jetbtagged", m_histjetbtagged) );

	m_histphoton = new TH1D("photon","photon",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/photon", m_histphoton) );
	m_histsignalphotonorg = new TH1D("signalphotonorg","signalphotonorg",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/signalphotonorg", m_histsignalphotonorg) );
	m_histsignalphoton = new TH1D("signalphoton","signalphoton",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/signalphoton", m_histsignalphoton) );
	m_histphotoncandidate = new TH1D("photoncandidate","photoncandidate",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/photoncandidate", m_histphotoncandidate) );
	m_histdiphoton = new TH1D("diphoton","diphoton",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/diphoton", m_histdiphoton) );

    /*for cutflow*/
	m_histpassedpt = new TH1D("passedpt","passedpt",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedpt", m_histpassedpt) );
	m_histpassedeta = new TH1D("passedeta","passedeta",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedeta", m_histpassedeta) );
	m_histpassedbadcls = new TH1D("passedbadcls","passedbadcls",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedbadcls", m_histpassedbadcls) );
	m_histpassedOQ = new TH1D("passedOQ","passedOQ",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedOQ", m_histpassedOQ) );
	m_histpassedauthor = new TH1D("passedauthor","passedauthor",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedauthor", m_histpassedauthor) );
	m_histpassedTightID = new TH1D("passedTightID","passedTightID",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedTightID", m_histpassedTightID) );
	m_histpassedTightID = new TH1D("passedTightID","passedTightID",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedTightID", m_histpassedTightID) );
	m_histpassedphoton = new TH1D("passedphoton","passedphoton",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/passedphoton", m_histpassedphoton) );
	m_histdiphotonbeforecut = new TH1D("diphotonbeforecut","diphotonbeforecut",500,0,500);
	CHECK( histSvc()->regHist("/MYSTREAM/diphotonbeforecut", m_histdiphotonbeforecut) );

    //tree
	m_tree = new TTree("tree","tree");
	CHECK( histSvc()->regTree("/MYSTREAM/tree", m_tree) );
	m_tree->Branch("runNumber", &m_runNumber);
	m_tree->Branch("eventNumber", &m_eventNumber);
	m_tree->Branch("lumiBlock", &m_lumiBlock);
	m_tree->Branch("averageIntPerXing", &m_averageIntPerXing);

	//jet
	TString jetAlgo = "AntiKt4EMTopo"; // Jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below) 
	TString config = "JES_MC16Recommendation_Consolidated_EMTopo_Apr2019_Rel21.config"; // Global config (see below) 
	TString calibSeq = "JetArea_Residual_EtaJES_GSC_Smear"; // Calibration sequence to apply (see below) 
	TString calibArea = "00-04-82"; // Calibration Area tag (see below) 
	bool isData = false; // bool describing if the events are data or from simulation 
	m_jetCalibration.setTypeAndName("JetCalibrationTool/JetCalibrationTool_EMTopo_MC");
	if( !m_jetCalibration.isUserConfigured() ){
		CHECK( m_jetCalibration.setProperty("JetCollection",jetAlgo.Data()) );
		CHECK( m_jetCalibration.setProperty("ConfigFile",config.Data()) );
		CHECK( m_jetCalibration.setProperty("CalibSequence",calibSeq.Data()) );
		CHECK( m_jetCalibration.setProperty("CalibArea",calibArea.Data()) );
		CHECK( m_jetCalibration.setProperty("IsData",isData) );
		CHECK( m_jetCalibration.retrieve() );
	}

	CHECK( m_btagSelectionTool.setProperty( "FlvTagCutDefinitionsFileName","xAODBTaggingEfficiency/13TeV/2017-21-13TeV-MC16-CDI-2018-02-09_v1.root" ) );
	CHECK( m_btagSelectionTool.setProperty("TaggerName", "MV2c10" ) );
	CHECK( m_btagSelectionTool.setProperty("OperatingPoint", "FixedCutBEff_70") );
	CHECK( m_btagSelectionTool.setProperty("JetAuthor", "AntiKt4EMTopoJets" ) );
	CHECK( m_btagSelectionTool.retrieve() );

	//photon
	CHECK( m_EgammaCalibrationAndSmearingTool.setProperty("ESModel", "es2018_R21_v0") );
	CHECK( m_EgammaCalibrationAndSmearingTool.retrieve() );
	CHECK( m_fudgeMCTool.setProperty("Preselection", 22) );
	CHECK( m_fudgeMCTool.setProperty("FFCalibFile","ElectronPhotonShowerShapeFudgeTool/v2/PhotonFudgeFactors.root") );
	CHECK( m_fudgeMCTool.retrieve() );
	//photon selection
	CHECK( m_photonLooseIsEMSelector.setProperty("WorkingPoint", "LoosePhoton") );
	CHECK( m_photonLooseIsEMSelector.retrieve() );
	CHECK( m_photonTightIsEMSelector.setProperty("WorkingPoint", "TightPhoton") );
	CHECK( m_photonTightIsEMSelector.retrieve() );
	CHECK( m_deadHVTool.retrieve() );
	//for diphoton
	CHECK( m_vertexTool.retrieve() );
	CHECK( m_pointTool.retrieve() );
	CHECK( m_indet_tracksel.setProperty("CutLevel" , "Loose") );
	CHECK( m_indet_tracksel.setProperty("minPt" , 1000.0 ) );
	CHECK( m_indet_tracksel.setProperty("maxZ0SinTheta", 3.0 ) );
	CHECK( m_indet_tracksel.retrieve() );
	CHECK( m_track_iso_tool.setProperty("TrackSelectionTool", m_indet_tracksel) );
	CHECK( m_track_iso_tool.retrieve() );
	bool isAtlfast = false;
	CHECK( m_isoCorrTool.setProperty( "IsMC", !isData) );
	CHECK( m_isoCorrTool.setProperty( "AFII_corr", false) );
	CHECK( m_isoCorrTool.retrieve() );
	CHECK( m_photonFCLooseIsoTool.setProperty("PhotonWP", "FixedCutLoose") );
	CHECK( m_photonFCLooseIsoTool.retrieve() );

	return StatusCode::SUCCESS;
}

StatusCode MyPackageAlg::finalize() {
	ATH_MSG_INFO ("Finalizing " << name() << "...");
	//
	//Things that happen once at the end of the event loop go here
	//


	return StatusCode::SUCCESS;
}

StatusCode MyPackageAlg::execute() {  
	ATH_MSG_DEBUG ("Executing " << name() << "...");
	setFilterPassed(false); //optional: start with algorithm not passed

	//Your main analysis code goes here
	//If you will use this algorithm to perform event skimming, you
	//should ensure the setFilterPassed method is called
	//If never called, the algorithm is assumed to have 'passed' by default
	// if data check if event passes GRL 
	bool pass_HLT_mu26_ivarmedium = m_tdt->isPassed("HLT_mu26_ivarmedium");
	bool pass_HLT_mu50 = m_tdt->isPassed("HLT_mu50");
	bool pass_HLT_g35loose = m_tdt->isPassed("HLT_g35_loose_g25_loose");
	bool pass_HLT_g35medium = m_tdt->isPassed("HLT_g35_medium_g25_medium_L12EM20VH");

	const xAOD::EventInfo* ei = 0;
	CHECK( evtStore()->retrieve( ei , "EventInfo" ) );
	//ATH_MSG_INFO("eventNumber=" << ei->eventNumber() );
	//m_myHist->Fill( ei->averageInteractionsPerCrossing() ); //fill mu into histogram

	// reject event if: 
	bool isMC = ei->eventType( xAOD::EventInfo::IS_SIMULATION );
	bool cleanEvent = true;

	if(!isMC){
		if( (ei->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) ||
				(ei->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) ||
				(ei->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) ||
				(ei->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) )
		{
			cleanEvent = false;
		} // end if event flags check 
	} // end if the event is data
	if(!cleanEvent) return StatusCode::SUCCESS;

	bool passGRL = true;
	if (!isMC) {
		if (!m_grl->passRunLB(*ei)) passGRL = false;
	} // end if the event is data
	if(!passGRL) return StatusCode::SUCCESS;

	m_runNumber = ei->runNumber();
	m_eventNumber = ei->eventNumber();
	m_lumiBlock = ei->lumiBlock();
	m_averageIntPerXing = ei->averageInteractionsPerCrossing();
	m_histAverageIntPerXing->Fill(m_averageIntPerXing);
	m_tree->Fill();

	//trigger check
	//if( !pass_HLT_mu26_ivarmedium || !pass_HLT_mu50 ) return StatusCode::SUCCESS;
	if( !pass_HLT_g35loose || !pass_HLT_g35medium ) return StatusCode::SUCCESS;

	//------------
	// MUONS
	//------------
	// get muon container of interest
	const xAOD::MuonContainer* muons = 0;
	CHECK (evtStore()->retrieve(muons, "Muons"));
	auto muons_shallowCopy = xAOD::shallowCopyContainer( *muons );
	std::unique_ptr<xAOD::MuonContainer> muonsSC (muons_shallowCopy.first);
	std::unique_ptr<xAOD::ShallowAuxContainer> muonsAuxSC (muons_shallowCopy.second);
	// loop over the muons in the container
	//for (auto muon : *muons) {
	//	  if ( !m_muonSelection->accept(*muon) ) continue;
	//	  ATH_MSG_INFO("execute(): original muon pt/eta/phi = " << ((muon)->pt() * 0.001) << " GeV/" << (muon)->eta() << "/" << (muon)->phi());
	// } // end for loop over muon
	// loop over the muons in the container 
	for (auto muonSC : *muonsSC) {
		muonSC->auxdata< int >( "passMedium" ) = m_muonSelection->accept(*muonSC);
		bool pass_iso = m_isoTool->accept(*muonSC);
		if(pass_iso){
			//ATH_MSG_INFO("execute(): original muon pt/eta/phi = " << (muonSC->pt() * 0.001) << " GeV/" << muonSC->eta() << "/" << muonSC->phi());
			if ( !m_muonSelection->accept(*muonSC) ) continue;
			if(m_muonCalibrationAndSmearingTool->applyCorrection(*muonSC) != CP::CorrectionCode::Ok){
				ATH_MSG_INFO ("execute(): Problem with Muon Calibration And Smearing Tool (Error or OutOfValidityRange) ");
			}
			//ATH_MSG_INFO("execute(): corrected muon pt/eta/phi = " << (muonSC->pt() * 0.001) << " GeV/" << muonSC->eta() << "/" << muonSC->phi());

			const xAOD::TrackParticle *IDtrack = muonSC->trackParticle(xAOD::Muon::TrackParticleType::InnerDetectorTrackParticle);
			if(IDtrack) {
				//ATH_MSG_INFO("execute(): ID track pt/eta/phi/d0 = " << IDtrack->pt() << "/" << IDtrack->eta() << "/" << IDtrack->phi() << "/" << IDtrack->d0());
			}
		}

		bool foundlink = false;
		ElementLink< xAOD::TruthParticleContainer > truthLink;
		if (muonSC->isAvailable< ElementLink< xAOD::TruthParticleContainer > > ("truthParticleLink") ){
			truthLink = muonSC->auxdata< ElementLink< xAOD::TruthParticleContainer> >("truthParticleLink");
			if(truthLink.isValid()){
				foundlink = true;
			}
			if (foundlink){
				//ATH_MSG_INFO("execute(): TruthParticle linked to reco muon barcode" << (*truthLink)->barcode() << " pdgId/status/pt = " << (*truthLink)->pdgId() << "/" << (*truthLink)->status() << "/" << (*truthLink)->pt());
			}
		}

	} // end for loop over muons

	//for dimuMass
	xAOD::MuonContainer::const_iterator mu_itr1 = muonsSC->begin();
	xAOD::MuonContainer::const_iterator mu_end = muonsSC->end();
	for ( ; mu_itr1!=mu_end; ++mu_itr1 ) {
		if(!(*mu_itr1)->auxdata< int >( "passMedium" )) continue;
		// second muon loop 
		xAOD::MuonContainer::const_iterator mu_itr2 = mu_itr1;
		for ( mu_itr2++; mu_itr2!=mu_end; ++mu_itr2 ) {
			if(!(*mu_itr2)->auxdata< int >( "passMedium" )) continue;
			if( (*mu_itr1)->charge()*(*mu_itr2)->charge() > 0 ) continue;
			// selection like if((*mu_itr1)->charge() * (*mu_itr2)->charge() > 0) continue;
			// then calculate di-muon mass
			const double dimuMass = ((*mu_itr1)->p4() + (*mu_itr2)->p4()).M();
			m_histdiMuon->Fill(dimuMass);
		}
	}

	//------------ 
	//// MC Truth 
	////------------ 
	if(isMC){
		const xAOD::TruthParticleContainer* mctruths = 0;
		CHECK ( evtStore()->retrieve( mctruths, "TruthParticles" ) );
		for (auto mctruth : *mctruths) {
			if(!mctruth) continue;
			float mctruth_pt = mctruth->pt();
			if(mctruth_pt<0.00001) continue;
			if(mctruth->pdgId() == 23 && mctruth->hasDecayVtx()) {
				//ATH_MSG_INFO("execute(): parent MC truth barcode" << mctruth->barcode() << " pdgId/status/pt = " << mctruth->pdgId() << "/" << mctruth->status() << "/" << mctruth_pt);
				//DumpDecayChain(mctruth);
			}
		}
	}


	//------------ 
	//// Jet 
	////------------ 
	// get jet container of interest 
	const xAOD::JetContainer* jets = 0;
	CHECK(evtStore()->retrieve( jets, "AntiKt4EMTopoJets" ));
	//ATH_MSG_INFO ("execute(): number of jets = " << jets->size());
	// loop over the jets in the container 
	//for (auto jet : *jets) {
	//	  ATH_MSG_INFO("execute(): original jet pt/eta/phi = " << (jet->pt() * 0.001) << " GeV/" << jet->eta() << "/" << jet->phi());
	// }
	// get jet container of interest 
	auto jets_shallowCopy = xAOD::shallowCopyContainer( *jets );
	std::unique_ptr<xAOD::JetContainer> jetsSC (jets_shallowCopy.first);
	std::unique_ptr<xAOD::ShallowAuxContainer> jetsAuxSC (jets_shallowCopy.second);
	// loop over the jets in the container 
	for (auto jetSC : *jetsSC) {
		// ATH_MSG_INFO ("execute(): jet pt = " << (jet->pt() * 0.001) << " GeV"); // just to print out something 
		//ATH_MSG_INFO("execute(): original jet pt/eta/phi = " << (jetSC->pt() * 0.001) << " GeV/" << jetSC->eta() << "/" << jetSC->phi());
		m_jetCalibration->applyCalibration(*jetSC);
		//ATH_MSG_INFO("execute(): corrected jet pt/eta/phi = " << (jetSC->pt() * 0.001) << " GeV/" << jetSC->eta() << "/" << jetSC->phi());

		bool tagged = m_btagSelectionTool->accept(*jetSC);
		double tagweight = 0;
		m_btagSelectionTool->getTaggerWeight( *jetSC ,tagweight);
		double cutvalue = 0;
		m_btagSelectionTool->getCutValue(jetSC->pt(), cutvalue);
		m_histjet->Fill(jetSC->pt()); 

		if(tagged) {
			//ATH_MSG_INFO("execute(): corrected jet pt/eta/phi = " << (jetSC->pt() * 0.001) << " GeV/" << jetSC->eta() << "/" << jetSC->phi());
			m_histjetbtagged->Fill(jetSC->pt()); 
		}

		if(isMC){
			int truthLabel = -1;
			jetSC->getAttribute("HadronConeExclExtendedTruthLabelID",truthLabel);
			if(tagged) ATH_MSG_INFO("execute(): truthlabelID = " << truthLabel);
		}

	} // end for loop over jets


	//------------ 
	// PHOTONS 
	//------------ 
	// get photon container of interest 
	const xAOD::PhotonContainer* photons = 0;
	CHECK (evtStore()->retrieve(photons, "Photons"));
	auto photons_shallowCopy = xAOD::shallowCopyContainer( *photons );
	std::unique_ptr<xAOD::PhotonContainer> photonsSC (photons_shallowCopy.first);
	std::unique_ptr<xAOD::ShallowAuxContainer> photonsAuxSC (photons_shallowCopy.second);

	bool PASS_PRESELECTION = false;

	// loop over the photons in the container 
	for (auto photonSC : *photonsSC) {
		m_histphoton->Fill(photonSC->pt() * 0.001); //pt(GeV)
		//ATH_MSG_INFO("execute(): original photon pt/eta/phi = " << (photonSC->pt() * 0.001) << " GeV/" << photonSC->eta() << "/" << photonSC->phi());
		if(m_EgammaCalibrationAndSmearingTool->applyCorrection(*photonSC) != CP::CorrectionCode::Ok){
			ATH_MSG_INFO ("execute(): Problem with Photon Calibration And Smearing Tool (Error or OutOfValidityRange)");
		}
		//ATH_MSG_INFO("execute(): corrected photon pt/eta/phi = " << (photonSC->pt() * 0.001) << " GeV/" << photonSC->eta() << "/" << photonSC->phi());

		if(isMC) { // NOT apply Fudge to AFII samples!!!!!!!!!!!!!!!!!!!!!!!!! 
			float unfudged_weta1 = photonSC->auxdata<float>("weta1");
			if(m_fudgeMCTool->applyCorrection(*photonSC) != CP::CorrectionCode::Ok){
				ATH_MSG_INFO ("execute(): Problem with Photon Shower Shift Correction (Error or OutOfValidityRange)");
			}
			float fudged_weta1 = photonSC->auxdata<float>("weta1");
			ATH_MSG_INFO ("execute(): weta1 = " << unfudged_weta1 << " -> " << fudged_weta1);
		}

		//select photon by eta
		bool pass_eta = false;
		if( fabs( photonSC->caloCluster()->etaBE(2) ) < 1.37
				|| (fabs( photonSC->caloCluster()->etaBE(2) ) > 1.52
					&& fabs( photonSC->caloCluster()->etaBE(2) ) < 2.37) ) {
			pass_eta = true;
		}

		//select photon by pt
		bool pass_pt = false;
		if (photonSC->pt()*0.001 > 25) pass_pt = true;

		//reject bad  photon 
		bool pass_badcls = photonSC->isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON);
		bool pass_OQ = false;
		if(!( (photonSC->OQ()&1073741824)!=0 ||
					( (photonSC->OQ()&134217728)!=0 &&
					  (photonSC->showerShapeValue(xAOD::EgammaParameters::Reta) >0.98
					   ||photonSC->showerShapeValue(xAOD::EgammaParameters::f1) > 0.4
					   ||(photonSC->OQ()&67108864)!=0)
					) ) ){
			pass_OQ = true;
		}

		//author cut
		bool pass_author = false;
		uint16_t author = photonSC->author();
		if ((author & xAOD::EgammaParameters::AuthorPhoton) || (author & xAOD::EgammaParameters::AuthorAmbiguous)){
			pass_author = true;
		}

		bool pass_LooseID = m_photonLooseIsEMSelector->accept(photonSC);
		bool pass_TightID = m_photonTightIsEMSelector->accept(photonSC);
		bool pass_HV = m_deadHVTool->accept(photonSC);

		ATH_MSG_INFO("execute(): photon pt/eta/phi = " << (photonSC->pt() * 0.001) << " GeV/" << photonSC->eta() << "/" << photonSC->phi());

		if( pass_pt && pass_eta && pass_badcls && pass_OQ && pass_author && pass_TightID /*&& pass_HV*/) PASS_PRESELECTION = true;

        //for cutflow
		if(pass_pt) m_histpassedpt->Fill(photonSC->pt() * 0.001); //pt(GeV)
		if(pass_eta) m_histpassedeta->Fill(photonSC->pt() * 0.001); //pt(GeV)
		if(pass_badcls) m_histpassedbadcls->Fill(photonSC->pt() * 0.001); //pt(GeV)
		if(pass_OQ) m_histpassedOQ->Fill(photonSC->pt() * 0.001); //pt(GeV)
		if(pass_author) m_histpassedauthor->Fill(photonSC->pt() * 0.001); //pt(GeV)
		if(pass_TightID) m_histpassedTightID->Fill(photonSC->pt() * 0.001); //pt(GeV)

		if(PASS_PRESELECTION) m_histpassedphoton->Fill(photonSC->pt() * 0.001); //pt(GeV)


	} // end for loop over photons

	const xAOD::Vertex *vertex = nullptr;
	xAOD::PhotonContainer signalphotons(photonsSC->begin(), photonsSC->end(), SG::VIEW_ELEMENTS);
	signalphotons.sort(comparePt);

	for (auto photon4vtx = signalphotons.begin(); photon4vtx != signalphotons.end();) {
		if( PASS_PRESELECTION )
			photon4vtx++;
		else
			photon4vtx = signalphotons.erase(photon4vtx);
	}

	if(signalphotons.size() < 2) {
		const xAOD::VertexContainer *vertices = nullptr;
		CHECK (evtStore()->retrieve(vertices, "PrimaryVertices"));
		vertex = xAOD::PVHelpers::getHardestVertex(vertices);
	}
	else {
		signalphotons.resize(2);
		m_vertexTool->getVertex(signalphotons, vertex);
	}

	for (auto signalphoton : signalphotons) {
		m_histsignalphotonorg->Fill(signalphoton->pt() * 0.001); //pt(GeV)

		if (signalphoton->caloCluster() && fabs(signalphoton->caloCluster()->etaBE(1)) < 10.0
				&& m_pointTool->correctPrimaryVertex(*signalphoton, vertex->z()).isFailure() ) {
			ATH_MSG_WARNING ("execute(): Pointing tool unable to correct photon for PV");
		}

		std::set<const xAOD::TrackParticle*> tracksToExclude;
		for (unsigned int iv = 0; iv < signalphoton->nVertices(); iv++) {
			const xAOD::Vertex *phvtx = signalphoton->vertex(iv);
			for (unsigned int itk = 0; itk < phvtx->nTrackParticles(); itk++) {
				tracksToExclude.insert(xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(phvtx->trackParticle(itk)));
			}
		}
		std::vector<xAOD::Iso::IsolationType> m_isoT = {xAOD::Iso::ptcone40, xAOD::Iso::ptcone30, xAOD::Iso::ptcone20};
		xAOD::TrackCorrection m_corrList;
		m_corrList.trackbitset.set(static_cast<unsigned int>(xAOD::Iso::coreTrackPtr));
		m_track_iso_tool->decorateParticle(*signalphoton, m_isoT, m_corrList, vertex, &tracksToExclude);
		if (not m_isoCorrTool->applyCorrection(*signalphoton))
			ATH_MSG_FATAL("execute(): Couldn't correct photon isolation leakage?");

		m_histsignalphoton->Fill(signalphoton->pt() * 0.001); //pt(GeV)

	}

	//for diphoton
	xAOD::PhotonContainer::const_iterator photon_itr1 = signalphotons.begin();
	xAOD::PhotonContainer::const_iterator photon_end = signalphotons.end();
	for ( ; photon_itr1!=photon_end; ++photon_itr1 ) {
		m_histphotoncandidate->Fill((*photon_itr1)->pt() * 0.001); //pt(GeV)
		// second photon loop 
		xAOD::PhotonContainer::const_iterator photon_itr2 = photon_itr1;
		for ( photon_itr2++; photon_itr2!=photon_end; ++photon_itr2 ) {
			// then get diphotonmass 
			const double diphotonmass = ((*photon_itr1)->p4() + (*photon_itr2)->p4()).M();
            if (diphotonmass*0.001 > 105 && diphotonmass*0.001 < 160) m_histdiphotonbeforecut->Fill(diphotonmass*0.001); //M(GeV)
            if ((*photon_itr1)->pt()/diphotonmass > 0.35 && (*photon_itr2)->pt()/diphotonmass > 0.25){
				if (diphotonmass*0.001 > 105 && diphotonmass*0.001 < 160) m_histdiphoton->Fill(diphotonmass*0.001); //M(GeV)
			}
		}
	}


	setFilterPassed(true); //if got here, assume that means algorithm passed
	return StatusCode::SUCCESS;
}

StatusCode MyPackageAlg::beginInputFile() { 
	//
	//This method is called at the start of each input file, even if
	//the input file contains no events. Accumulate metadata information here
	//

	//example of retrieval of CutBookkeepers: (remember you will need to include the necessary header files and use statements in requirements file)
	// const xAOD::CutBookkeeperContainer* bks = 0;
	// CHECK( inputMetaStore()->retrieve(bks, "CutBookkeepers") );

	//example of IOVMetaData retrieval (see https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AthAnalysisBase#How_to_access_file_metadata_in_C)
	//float beamEnergy(0); CHECK( retrieveMetadata("/TagInfo","beam_energy",beamEnergy) );
	//std::vector<float> bunchPattern; CHECK( retrieveMetadata("/Digitiation/Parameters","BeamIntensityPattern",bunchPattern) );



	return StatusCode::SUCCESS;
}

StatusCode MyPackageAlg::DumpDecayChain(const xAOD::TruthParticle* parent) {
	bool hasDecayVtx = parent->hasDecayVtx();
	if(!hasDecayVtx) return StatusCode::SUCCESS;
	const xAOD::TruthVertex* vtx = parent->decayVtx();
	float vtx_x = vtx->x();
	float vtx_y = vtx->y();
	float vtx_z = vtx->z();
	float vtx_r = sqrt(vtx_x*vtx_x + vtx_y*vtx_y);
	float vtx_phi = atan2(vtx_x,vtx_y);
	ATH_MSG_INFO("DumpDecayChain(): barcode" << parent->barcode() << " decays at Z/R/phi = " << vtx_z << "/" << vtx_r << "/" << vtx_phi);
	for(unsigned int i=0; i<vtx->nOutgoingParticles(); i++) {
		const xAOD::TruthParticle* child = vtx->outgoingParticle(i);
		if(!child) continue;
		float child_pt = child->pt();
		if(child_pt<0.00001) continue;
		ATH_MSG_INFO("DumpDecayChain(): child of " << parent->barcode() << ", barcode" << child->barcode() << " pdgId/status/pt = " << child->pdgId() << "/" << child->status() << "/" << child_pt);
		DumpDecayChain(child);
	}
	return StatusCode::SUCCESS;
}

bool MyPackageAlg::comparePt(const xAOD::Photon *a, const xAOD::Photon *b)
{
	return a->pt() > b->pt();
}

