#ifndef MYPACKAGE_MYPACKAGEALG_H
#define MYPACKAGE_MYPACKAGEALG_H 1

#include "AthAnalysisBaseComps/AthAnalysisAlgorithm.h"
#include "AsgTools/AnaToolHandle.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h"
#include <MuonAnalysisInterfaces/IMuonSelectionTool.h>
#include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
#include "IsolationSelection/IIsolationSelectionTool.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "JetCalibTools/IJetCalibrationTool.h"
#include "xAODBTaggingEfficiency/BTaggingSelectionTool.h"
#include "EgammaAnalysisInterfaces/IEgammaCalibrationAndSmearingTool.h"
#include "EgammaAnalysisInterfaces/IElectronPhotonShowerShapeFudgeTool.h"
#include "EgammaAnalysisInterfaces/IAsgPhotonIsEMSelector.h"
#include "EgammaAnalysisInterfaces/IAsgDeadHVCellRemovalTool.h"
#include "PhotonVertexSelection/IPhotonVertexSelectionTool.h"
#include "PhotonVertexSelection/IPhotonPointingTool.h"
#include "RecoToolInterfaces/ITrackIsolationTool.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "IsolationCorrections/IIsolationCorrectionTool.h"


//Example ROOT Includes
//#include "TTree.h"
//#include "TH1D.h"



class MyPackageAlg: public ::AthAnalysisAlgorithm { 
 public: 
  MyPackageAlg( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~MyPackageAlg(); 

  ///uncomment and implement methods as required

                                        //IS EXECUTED:
  virtual StatusCode  initialize();     //once, before any input is loaded
  virtual StatusCode  beginInputFile(); //start of each input file, only metadata loaded
  //virtual StatusCode  firstExecute();   //once, after first eventdata is loaded (not per file)
  virtual StatusCode  execute();        //per event
  //virtual StatusCode  endInputFile();   //end of each input file
  //virtual StatusCode  metaDataStop();   //when outputMetaStore is populated by MetaDataTools
  virtual StatusCode  DumpDecayChain(const xAOD::TruthParticle* parent);
  virtual StatusCode  finalize();       //once, after all events processed
  //for di-photon
  static bool comparePt(const xAOD::Photon *a, const xAOD::Photon *b);
  

  ///Other useful methods provided by base class are:
  ///evtStore()        : ServiceHandle to main event data storegate
  ///inputMetaStore()  : ServiceHandle to input metadata storegate
  ///outputMetaStore() : ServiceHandle to output metadata storegate
  ///histSvc()         : ServiceHandle to output ROOT service (writing TObjects)
  ///currentFile()     : TFile* to the currently open input file
  ///retrieveMetadata(...): See twiki.cern.ch/twiki/bin/view/AtlasProtected/AthAnalysisBase#ReadingMetaDataInCpp


 private: 

   //Example algorithm property, see constructor for declaration:
   //int m_nProperty = 0;

   //Example histogram, see initialize method for registration to output histSvc
   //TH1D* m_myHist = 0;
   //TTree* m_myTree = 0;

  //histograms
  TH1D* m_histAverageIntPerXing = 0;
  TH1D* m_histdiMuon = 0;
  TH1D* m_histjet = 0;
  TH1D* m_histjetbtagged = 0;
  TH1D* m_histphoton = 0;
  TH1D* m_histpassedphoton = 0;
  TH1D* m_histsignalphotonorg = 0;
  TH1D* m_histsignalphoton = 0;
  TH1D* m_histphotoncandidate = 0;
  TH1D* m_histdiphoton = 0;
  /*for cutflow*/
  m_histpassedpt = 0; 
  m_histpassedeta = 0; 
  m_histpassedbadcls = 0; 
  m_histpassedOQ = 0; 
  m_histpassedauthor = 0; 
  m_histpassedTightID = 0;
  m_histpassedTightID = 0;
  m_histpassedphoton = 0; 
  m_histdiphotonbeforecut = 0;

  //tree
  TTree* m_tree = 0;

  uint32_t m_runNumber = 0;
  unsigned long long m_eventNumber = 0;
  uint32_t m_lumiBlock = 0;
  float m_averageIntPerXing = 0;
  
  asg::AnaToolHandle<Trig::TrigDecisionTool> m_tdt;
  asg::AnaToolHandle<IGoodRunsListSelectionTool> m_grl;
  asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelection;
  asg::AnaToolHandle<CP::IMuonCalibrationAndSmearingTool> m_muonCalibrationAndSmearingTool;
  asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoTool;
  asg::AnaToolHandle<IJetCalibrationTool> m_jetCalibration;
  asg::AnaToolHandle<IBTaggingSelectionTool> m_btagSelectionTool;
  asg::AnaToolHandle<CP::IEgammaCalibrationAndSmearingTool> m_EgammaCalibrationAndSmearingTool;
  asg::AnaToolHandle<IElectronPhotonShowerShapeFudgeTool> m_fudgeMCTool;
  asg::AnaToolHandle<IAsgPhotonIsEMSelector> m_photonLooseIsEMSelector;
  asg::AnaToolHandle<IAsgPhotonIsEMSelector> m_photonTightIsEMSelector;
  asg::AnaToolHandle<IAsgDeadHVCellRemovalTool> m_deadHVTool;
  //for di-photon
  asg::AnaToolHandle<CP::IPhotonVertexSelectionTool> m_vertexTool;
  asg::AnaToolHandle<CP::IPhotonPointingTool> m_pointTool;
  asg::AnaToolHandle<InDet::IInDetTrackSelectionTool> m_indet_tracksel;
  asg::AnaToolHandle<xAOD::ITrackIsolationTool> m_track_iso_tool;
  asg::AnaToolHandle<CP::IIsolationCorrectionTool> m_isoCorrTool;
  asg::AnaToolHandle<CP::IIsolationSelectionTool> m_photonFCLooseIsoTool;

}; 

#endif //> !MYPACKAGE_MYPACKAGEALG_H
