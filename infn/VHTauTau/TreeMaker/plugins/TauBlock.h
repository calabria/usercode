#ifndef __VHTauTau_TreeMaker_TauBlock_h
#define __VHTauTau_TreeMaker_TauBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Ref.h"

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"

class TClonesArray;
class Tau;

class TauBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TauBlock(const edm::ParameterSet& iConfig);
  virtual ~TauBlock();

  static const reco::PFJetRef& getJetRef(const reco::PFTau& tau);

  enum {
    kMaxTau = 100
  };

private:
  TClonesArray* cloneTau; 
  int  fnTau;
  int _verbosity;
  edm::InputTag _inputTag;
  edm::InputTag _vtxInputTag;
  edm::InputTag muonSrc_;
  edm::InputTag srcBeamSpot_;
  edm::InputTag _beamSpotInputTag;
  bool _beamSpotCorr;
  edm::InputTag genParticleSrc_;
  edm::ParameterSet qualityCutsPSet_;
  edm::InputTag pfCandSrc_;
  edm::InputTag vertexSrc_;

  std::auto_ptr<reco::tau::RecoTauVertexAssociator> vertexAssociator_;
  std::auto_ptr<reco::tau::RecoTauQualityCuts> qcuts_;
  // Inverted QCut which selects tracks with bad DZ/trackWeight
  std::auto_ptr<reco::tau::RecoTauQualityCuts> pileupQcutsPUTrackSelection_;
  std::auto_ptr<reco::tau::RecoTauQualityCuts> pileupQcutsGeneralQCuts_;

  typedef std::vector<double> vdouble;

  const TransientTrackBuilder* trackBuilder_;

  vhtm::Tau* tauB;
};
#endif
