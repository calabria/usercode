import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
  verbosity       = cms.int32(0),
  patTauSrc       = cms.InputTag('tauVariables'),
  vertexSrc       = cms.InputTag('goodOfflinePrimaryVertices'),
  beamSpotCorr    = cms.bool(True),
  muonSrc         = cms.InputTag("selectedMuonsPt"),
  srcBeamSpot     = cms.untracked.InputTag("offlineBeamSpot"),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  genParticles    = cms.InputTag("genParticles"),
  qualityCuts = PFTauQualityCuts,# set the standard quality cuts
  #Delta-Beta corrections to remove Pileup
  particleFlowSrc = cms.InputTag("particleFlow"),
)
