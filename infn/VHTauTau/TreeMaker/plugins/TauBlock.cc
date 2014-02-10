#include <iostream>
#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Utilities/General/interface/FileInPath.h"

#include "VHTauTau/TreeMaker/plugins/TauBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/Ref.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <functional>
#include <boost/foreach.hpp>
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/ConeTools.h"

using namespace reco;
using namespace std;
using namespace edm;

namespace tb 
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref)
  {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}

//-------------------------------------------------------------------------------
 // auxiliary functions to compare two Tracks 


bool tracksMatchByDeltaR(const reco::Track* trk1, const reco::Track* trk2)
{
  if ( reco::deltaR(*trk1, *trk2) < 1.e-2 && trk1->charge() == trk2->charge() ) return true;
  else return false;
}

// auxiliary function to exclude tracks associated to tau lepton decay "leg"
// from primary event vertex refit
typedef std::map<const reco::Track*, reco::TransientTrack> TransientTrackMap;
void removeTracks(TransientTrackMap& pvTracks_toRefit, const std::vector<reco::Track*> svTracks)
{
  for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin(); svTrack != svTracks.end(); ++svTrack ) {
    
    //--- remove track from list of tracks included in primary event vertex refit
    //    if track matches by reference or in eta-phi
    //    any of the tracks associated to tau lepton decay "leg"
    for ( TransientTrackMap::iterator pvTrack = pvTracks_toRefit.begin(); pvTrack != pvTracks_toRefit.end(); ++pvTrack ) {
      if ( tracksMatchByDeltaR(pvTrack->first, *svTrack) ) {
	pvTracks_toRefit.erase(pvTrack);
	break;
      }
    }
  }
}

//-------------------------------------------------------------------------------

double vertexSignificance(reco::Vertex& pv, reco::Vertex& sv, GlobalVector& direction){

  return SecondaryVertex::computeDist3d(pv,sv,direction,true).significance();

}

//-------------------------------------------------------------------------------

std::vector<double> GenSecondaryVertex(std::vector<pat::Tau>::const_iterator it, edm::Handle<reco::GenParticleCollection> genParticles){

  std::vector<double> vtxCoord;

  //int runNumber = iEvent.id().run();
  double phiReco = it->phi();
  double etaReco = it->eta();

  if(genParticles.isValid()){

  	for(reco::GenParticleCollection::const_iterator itg = genParticles->begin(); itg != genParticles->end(); ++itg ){

		int id = itg->pdgId();
		int status = itg->status();
		//std::cout<<"Id = "<<id<<std::endl;
		int nDaughters = itg->numberOfDaughters();
		double phiGen = itg->phi();
		double etaGen = itg->eta();

		//std::cout<<"id "<<id<<" "<<phiGen<<" "<<etaGen<<std::endl;

		if(fabs(id) == 23 && status == 3){

			vtxCoord.push_back(itg->vx()); 
			vtxCoord.push_back(itg->vy());
			vtxCoord.push_back(itg->vz());

		}

		int numNeutrinoTau = 0;
		int numNeutrinos = 0;

		if(abs(id) == 15){

			for ( int k = 0; k < nDaughters; k++ ){

				const reco::Candidate * DaughterTau = itg->daughter(k);

				int idDaughterTau = DaughterTau->pdgId();
				//std::cout<<"idDaughterTau = "<<idDaughterTau<<std::endl;
				if(fabs(idDaughterTau) == 16) ++numNeutrinoTau;
				if(fabs(idDaughterTau) == 12 || fabs(idDaughterTau) == 14) ++numNeutrinos;

			}

		}
		else{

			numNeutrinoTau = 1;
			numNeutrinos = 0;

		}

		if(abs(id) == 15 && numNeutrinoTau == 1 && numNeutrinos == 0){

			//std::cout<<"Tutto "<<phiGen<<" "<<etaGen<<" "<<etaReco<<" "<<phiReco<<std::endl;

			double dR = deltaR(etaGen, phiGen, etaReco, phiReco);

			if(dR <= 0.5 && numNeutrinoTau == 1 && numNeutrinos == 0){

				const reco::Candidate * DaughterTau = itg->daughter(0);

				vtxCoord.push_back(DaughterTau->vx()); 
				vtxCoord.push_back(DaughterTau->vy());
				vtxCoord.push_back(DaughterTau->vz());
				vtxCoord.push_back(itg->px()); 
				vtxCoord.push_back(itg->py());
				vtxCoord.push_back(itg->pz());
				//std::cout<<"Id = "<<id<<" dR = "<<dR<<" numNeutrinoTau = "<<numNeutrinoTau<<std::endl;
				//std::cout<<" numNeutrinos = "<<numNeutrinos<<std::endl;

			}

		}

	}

  }

  return vtxCoord;

}

//-------------------------------------------------------------------------------

std::vector<double> produceNewVars(std::vector<pat::Tau>::const_iterator it, edm::Handle<edm::View<pat::Muon> > muons, edm::Handle<reco::VertexCollection> pvHandle, edm::Handle<reco::BeamSpot> beamSpotHandle, edm::ESHandle<TransientTrackBuilder> trackBuilderHandle){

   std::vector<double> newVars;
   VertexState beamSpotState_;
   bool beamSpotIsValid_;
   const reco::BeamSpot*  beamSpot_ = beamSpotHandle.product();

   if ( beamSpot_ ) {
     beamSpotState_ = VertexState(*beamSpot_);    
     beamSpotIsValid_ = (beamSpotState_.error().cxx() > 0. &&
			 beamSpotState_.error().cyy() > 0. &&
			 beamSpotState_.error().czz() > 0.);
   } else {
     edm::LogError ("NSVfitEventVertexRefitter::beginEvent")
       << " Failed to access BeamSpot !!";
     beamSpotIsValid_ = false;
   }
   
   const TransientTrackBuilder* trackBuilder_;
   trackBuilder_ = trackBuilderHandle.product();

   std::vector<reco::TransientTrack> tracks;

   //if (it->pfJetRef().isAvailable() && it->pfJetRef().isNonnull()) std::cout<<"pfJet: "<<it->pfJetRef()->pt()<<std::endl;
   //if (it->pfJetRef().isAvailable() && it->pfJetRef().isNonnull()) std::cout<<"pfJet: "<<it->pfJetRef()->phi()<<std::endl;
   //if (it->pfJetRef().isAvailable() && it->pfJetRef().isNonnull()) std::cout<<"pfJet: "<<it->pfJetRef()->eta()<<std::endl;

   const reco::PFCandidateRefVector& signalPFChargedHadrCandidates = it->signalPFChargedHadrCands();

   double SecVtx_x = -9, SecVtx_y = -9, SecVtx_z = -9;
   double PV_x = -9, PV_y = -9, PV_z = -9;
   //bool isSV = false;

   int isOneProng = 9999;
   int isThreeProng = 9999;
   double IPxyLead = 9999;
   double sigma_IPxyLead = 9999;
   double IPxySigLead = 9999;
   double IPzLead = 9999;

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////// Zona IP

   if((it->signalPFGammaCands().size() > 0 || it->signalPFGammaCands().size() == 0) && it->signalPFChargedHadrCands().size() == 1){

 	isOneProng = 1;
	isThreeProng = 0;

   }
   else if(signalPFChargedHadrCandidates.size() == 3){

 	isOneProng = 0;
	isThreeProng = 1;

   }
   else{

 	isOneProng = 0;
	isThreeProng = 0;
   }

   std::vector<reco::TransientTrack> pvTracks_original;
   TransientTrackMap pvTrackMap_refit;
       
   for ( reco::Vertex::trackRef_iterator pvTrack = (*pvHandle)[0].tracks_begin(); pvTrack != (*pvHandle)[0].tracks_end(); ++pvTrack ) {
	 
	reco::TransientTrack pvTrack_transient = trackBuilder_->build(pvTrack->get());
	pvTracks_original.push_back(pvTrack_transient);
	pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient)); //tracce associate al PV
	 
   }

   if (isOneProng){
     
     	std::vector<reco::Track*> svTracks;
     
     	if(it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()){
       
       		reco::Track trk_tmp0 = *(it->signalPFChargedHadrCands()[0]->trackRef());
       		reco::Track* ftrk_tmp0 = &trk_tmp0; 
       		svTracks.push_back(ftrk_tmp0);

	}
       
     
       	for(edm::View<pat::Muon>::const_iterator patMu=muons->begin(); patMu!=muons->end(); ++patMu){
	 
	 	if(patMu->track().isNonnull()){
	 
	   		reco::Track trk_tmp_mu = *(patMu->track());
	   		reco::Track* ftrk_tmp_mu = &trk_tmp_mu; 
	   		svTracks.push_back(ftrk_tmp_mu);

		}
       	}

   	removeTracks(pvTrackMap_refit, svTracks); //rimuovo dalla mappa gli elementi in cui le tracce associate al PV coincidono con quelle associate al SV

	    ////////////////definisco un vettore con le tracce del PV da cui ho tolto le tracce del SV///////////////////////////////////////////

   	std::vector<reco::TransientTrack> pvTracks_refit;
   	for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {

		pvTracks_refit.push_back(pvTrack->second);

   	}


	AdaptiveVertexFitter adf;
	TransientVertex pvVtx;

	if ( pvTracks_refit.size() >= 2) {
	      
		if (  beamSpotIsValid_ )  pvVtx = adf.vertex(pvTracks_refit, *beamSpot_);
		else pvVtx = adf.vertex(pvTracks_refit);

	}
	else{

		if ( beamSpotIsValid_ ) pvVtx = adf.vertex(pvTracks_original, *beamSpot_);
		else pvVtx = adf.vertex(pvTracks_original);

	}

	GlobalPoint PV =  pvVtx.position();
	PV_x = PV.x();
	PV_y = PV.y();
	PV_z = PV.z();
	/*std::cout<<"sono nel treeMaker, 1prong PVRef PV_X="<< PV_x<<std::endl; 
	std::cout<<"sono nel treeMaker, 1prong PVRef PV_Y="<< PV_y<<std::endl; 
	std::cout<<"sono nel treeMaker, 1prong PVRef PV_Z="<< PV_z<<std::endl;*/

	reco::Vertex::Point pvPoint = reco::Vertex::Point(pvVtx.position());
	    
	//std::cout<<"sono nel treeMaker, impact parameter wrt PVREFIT "<<it->leadPFChargedHadrCand()->trackRef()->dxy(pvPoint)<<"  err    "<<it->leadPFChargedHadrCand()->trackRef()->d0Error()<<std::endl;

	IPxyLead = it->leadPFChargedHadrCand()->trackRef()->dxy(pvPoint);
	sigma_IPxyLead = it->leadPFChargedHadrCand()->trackRef()->d0Error();
	IPxySigLead = IPxyLead/sigma_IPxyLead;
	IPzLead = it->leadPFChargedHadrCand()->trackRef()->dz(pvPoint);

   }
   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////// Zona SV

   reco::TransientTrack transtrack0; 
   reco::TransientTrack transtrack1; 
   reco::TransientTrack transtrack2; 

   double distSV_PVRef = -999;
   double flightPathSignificance = -999;
    
   if(isThreeProng){

	  if(it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()) transtrack0 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[0]->trackRef());
	  if(it->signalPFChargedHadrCands()[1]->trackRef().isNonnull()) transtrack1 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[1]->trackRef());
	  if(it->signalPFChargedHadrCands()[2]->trackRef().isNonnull()) transtrack2 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[2]->trackRef());
	  if (transtrack0.isValid()) tracks.push_back(transtrack0);
	  if (transtrack1.isValid()) tracks.push_back(transtrack1);
	  if (transtrack2.isValid()) tracks.push_back(transtrack2);
	  
	  //float vtxChi2 = -1;
	  //float vtxNDOF = -1;

	  if (tracks.size() > 1) {

	    	KalmanVertexFitter kvf(true);
	    	TransientVertex vtx = kvf.vertex(tracks);
	    	//vtxChi2 = vtx.totalChiSquared();
	    	//vtxNDOF = vtx.degreesOfFreedom();
	    
	    	//std::cout<<" vtx chi2="<<vtxChi2<<std::endl;
	    	//std::cout<<" vtx NDOF="<<vtxNDOF<<std::endl;
	    	//std::cout<<" vtx pos="<<vtx.position()<<std::endl;

	    	GlobalPoint SecVtx = vtx.position();
	    	SecVtx_x =  SecVtx.x();
	    	SecVtx_y =  SecVtx.y();
	    	SecVtx_z =  SecVtx.z();

		/*std::cout<<"3 prong SV_X="<< SecVtx_x<<std::endl; 
		std::cout<<"3 prong SV_Y="<< SecVtx_y<<std::endl; 
		std::cout<<"3 prong SV_Z="<< SecVtx_z<<std::endl;*/
	 
	    	///////////////////////////metto le tracce del 3 prong + la traccia del mu in un vettore//////////////////////////////////////////////////

	    	std::vector<reco::Track*> svTracks;
	    	bool isValidTauTrack = false;
	    	if((it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()) && (it->signalPFChargedHadrCands()[1]->trackRef().isNonnull()) && (it->signalPFChargedHadrCands()[2]->trackRef().isNonnull()) ) isValidTauTrack=true;

	    	if(isValidTauTrack){

	    		reco::Track trk_tmp0 = *(it->signalPFChargedHadrCands()[0]->trackRef());
	    		reco::Track trk_tmp1 = *(it->signalPFChargedHadrCands()[1]->trackRef());
	    		reco::Track trk_tmp2 = *(it->signalPFChargedHadrCands()[2]->trackRef());
	    		reco::Track* ftrk_tmp0 = &trk_tmp0; 
	    		reco::Track* ftrk_tmp1 = &trk_tmp1; 
	    		reco::Track* ftrk_tmp2 = &trk_tmp2; 
	    		svTracks.push_back(ftrk_tmp0);
	    		svTracks.push_back(ftrk_tmp1);
	    		svTracks.push_back(ftrk_tmp2);

	    	}

	    	for(edm::View<pat::Muon>::const_iterator patMu=muons->begin(); patMu!=muons->end(); ++patMu){

			if(patMu->track().isNonnull()){

		   		reco::Track trk_tmp_mu = *(patMu->track());
		   		reco::Track* ftrk_tmp_mu = &trk_tmp_mu; 
		   		svTracks.push_back(ftrk_tmp_mu);

			}
		 
           	}
	    
	    	//std::cout<<" puntatore  traccia "<<ftrk_tmp0<<"  =  "<<&trk_tmp0<<"  pt   "<<(*svTracks[0]).pt()<<std::endl;

	    	/*int idx=0;
	     	for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin();  svTrack != svTracks.end(); ++svTrack ) {
	      		std::cout << "Track #" << idx << ": Pt = " << (*svTrack)->pt() << "," 
			<< " eta = " << (*svTrack)->eta() << ", phi = " << (*svTrack)->phi() 
			<< " (charge = " << (*svTrack)->charge() << ", chi2 = " << (*svTrack)->normalizedChi2() << ")" << std::endl;
	      	      	++idx;
		}*/

	    	//isSV = true;
	    	removeTracks(pvTrackMap_refit, svTracks); //rimuovo dalla mappa gli elementi in cui le tracce associate al PV coincidono con quelle associate al SV

	    	////////////////definisco un vettore con le tracce del PV da cui ho tolto le tracce del SV///////////////////////////////////////////

	    	std::vector<reco::TransientTrack> pvTracks_refit;
	    	for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {

	      		pvTracks_refit.push_back(pvTrack->second);

	    	}

	    	//////////////////////////////////////////Refit del PV//////////////////////////////////////////////////////////////////////////////

	    	AdaptiveVertexFitter adf;
	    	TransientVertex pvVtx;

	    	if ( pvTracks_refit.size() >= 2) {
	      
	      		if (  beamSpotIsValid_ )  pvVtx = adf.vertex(pvTracks_refit, *beamSpot_);
	      		else pvVtx = adf.vertex(pvTracks_refit);

	    	}
	    	else{

	      		if ( beamSpotIsValid_ ) pvVtx = adf.vertex(pvTracks_original, *beamSpot_);
	      		else pvVtx = adf.vertex(pvTracks_original);

	    	}

	    	GlobalPoint PV =  pvVtx.position();
	    	PV_x = PV.x();
	    	PV_y = PV.y();
		PV_z = PV.z();
		/*std::cout<<"3 prong PV_X="<< PV_x<<std::endl; 
		std::cout<<"3 prong PV_Y="<< PV_y<<std::endl; 
		std::cout<<"3 prong PV_Z="<< PV_z<<std::endl;*/

    		if(vtx.isValid()){

      			GlobalVector tauDir(it->px(), it->py(), it->pz());
      			reco::Vertex secVer = vtx;
      			reco::Vertex primaryVertexNonConst = pvVtx;
      			flightPathSignificance = vertexSignificance(primaryVertexNonConst,secVer,tauDir);
			std::cout<<"flightPathSignificance "<<flightPathSignificance<<std::endl;

    		}

		distSV_PVRef = TMath::Sqrt((PV_x-SecVtx_x)*(PV_x-SecVtx_x) + (PV_y-SecVtx_y)*(PV_y-SecVtx_y) + (PV_z-SecVtx_z)*(PV_z-SecVtx_z)); 
		//std::cout<<"  ***dist pvRefit-sv***  "<<distSV_PVRef<<std::endl;

		reco::Vertex::Point pvPoint = reco::Vertex::Point(pvVtx.position());
	    
		//std::cout<<"sono nel treeMaker, impact parameter wrt PVREFIT "<<it->leadPFChargedHadrCand()->trackRef()->dxy(pvPoint)<<"  err    "<<it->leadPFChargedHadrCand()->trackRef()->d0Error()<<std::endl;

		IPxyLead = it->leadPFChargedHadrCand()->trackRef()->dxy(pvPoint);
		sigma_IPxyLead = it->leadPFChargedHadrCand()->trackRef()->d0Error();
		IPxySigLead = IPxyLead/sigma_IPxyLead;
		IPzLead = it->leadPFChargedHadrCand()->trackRef()->dz(pvPoint);

	  }

   }//is 3 prong

   newVars.push_back(isOneProng);
   newVars.push_back(isThreeProng);
   newVars.push_back(IPxyLead);
   newVars.push_back(sigma_IPxyLead);
   newVars.push_back(IPxySigLead);
   newVars.push_back(IPzLead);
   newVars.push_back(distSV_PVRef);
   newVars.push_back(PV_x);
   newVars.push_back(PV_y);
   newVars.push_back(PV_z);
   newVars.push_back(SecVtx_x);
   newVars.push_back(SecVtx_y);
   newVars.push_back(SecVtx_z);
   newVars.push_back(flightPathSignificance);
	
   /*if(isSV){

	//std::cout<<" sec vtx posX="<<SecVtx_x<<std::endl;
	if (!pvHandle->empty()) {

		const reco::Vertex &pv = (*pvHandle)[0];
		float pvx = pv.x();
		float pvy = pv.y();
		float pvz = pv.z();
		//std::cout<<" PV_X="<<pv.x()<<std::endl; 
		// std::cout<<" PV_Y="<<pv.y()<<std::endl; 
		//std::cout<<" PV_Z="<<pv.z()<<std::endl; 
		float distPS = TMath::Sqrt((pvx-SecVtx_x)*(pvx-SecVtx_x) + (pvy-SecVtx_y)*(pvy-SecVtx_y) + (pvz-SecVtx_z)*(pvz-SecVtx_z));
		//std::cout<<"  ***dist pv-sv***  "<< distPS<<std::endl;

	}

   }*/
	 
   return newVars;

}

//-------------------------------------------------------------------------------

std::vector<double> produceIso(std::vector<pat::Tau>::const_iterator it, edm::Handle<reco::PFCandidateCollection> pfCandHandle_, std::auto_ptr<tau::RecoTauQualityCuts> qcuts_, std::auto_ptr<tau::RecoTauQualityCuts> pileupQcutsPUTrackSelection_, std::auto_ptr<tau::RecoTauQualityCuts> pileupQcutsGeneralQCuts_, double coneRadius = 0.5, double dBCone = 0.8, double deltaBetaFactor = 0.3){

  std::vector<double> results;

  std::vector<reco::PFCandidateRef> chargedPFCandidatesInEvent_;
  std::vector<reco::PFCandidateRef> chargedHadPFCandidatesInEvent_;
  std::vector<reco::PFCandidateRef> gammaPFCandidatesInEvent_;
  std::vector<reco::PFCandidateRef> isoPU;

  for (size_t i = 0; i < pfCandHandle_->size(); ++i) {

      	reco::PFCandidateRef pfCand(pfCandHandle_, i);
	if (pfCand->charge() != 0) chargedPFCandidatesInEvent_.push_back(pfCand);
      	if (PFCandidate::ParticleType((*pfCand).particleId()) == PFCandidate::h) chargedHadPFCandidatesInEvent_.push_back(pfCand);
	else if(PFCandidate::ParticleType((*pfCand).particleId()) == PFCandidate::gamma) gammaPFCandidatesInEvent_.push_back(pfCand);

  }

  typedef reco::tau::cone::DeltaRPtrFilter<PFCandidateRef> DRFilter;
  DRFilter filter(it->p4(), 0, coneRadius);

  std::vector<PFCandidateRef> isoCharged_filter;
  std::vector<PFCandidateRef> isoNeutral_filter;

  BOOST_FOREACH(const PFCandidateRef& isoObject, chargedHadPFCandidatesInEvent_) {

      if(filter(isoObject)) isoCharged_filter.push_back(isoObject);

  }
  BOOST_FOREACH(const PFCandidateRef& isoObject, gammaPFCandidatesInEvent_) {

      if(filter(isoObject)) isoNeutral_filter.push_back(isoObject);

  }

  // First select by inverted the DZ/track weight cuts. True = invert
  //std::cout << "Initial PFCands: " << chargedPFCandidatesInEvent_.size()
  //  << std::endl;

  std::vector<PFCandidateRef> allPU = pileupQcutsPUTrackSelection_->filterRefs(chargedPFCandidatesInEvent_, true);

  // Now apply the rest of the cuts, like pt, and TIP, tracker hits, etc
  std::vector<PFCandidateRef> cleanPU = pileupQcutsGeneralQCuts_->filterRefs(allPU);

  // Only select PU tracks inside the isolation cone.
  DRFilter deltaBetaFilter(it->p4(), 0, dBCone);
  BOOST_FOREACH(const reco::PFCandidateRef& cand, cleanPU) {

      if (deltaBetaFilter(cand)) isoPU.push_back(cand);

  }

  double chargedPt = 0.0;
  double puPt = 0.0;
  double neutralPt = 0.0;

  BOOST_FOREACH(const PFCandidateRef& isoObject, isoCharged_filter) {

      	chargedPt += isoObject->pt();

  }
  BOOST_FOREACH(const PFCandidateRef& isoObject, isoNeutral_filter) {

      	neutralPt += isoObject->pt();

  }
  BOOST_FOREACH(const PFCandidateRef& isoObject, isoPU) {

      	puPt += isoObject->pt();

  }

  results.push_back(chargedPt);
  results.push_back(neutralPt);
  results.push_back(puPt);

  std::cout<<"Ch: "<<results[0]<<" Gamma: "<<results[1]<<" dB: "<<results[2]<<std::endl;

  neutralPt -= deltaBetaFactor*puPt;

  if (neutralPt < 0.0) neutralPt = 0.0;

  results.push_back(neutralPt);

  return results;

}

//-------------------------------------------------------------------------------

TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  muonSrc_(iConfig.getParameter<edm::InputTag>("muonSrc")),
  _beamSpotInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _beamSpotCorr(iConfig.getParameter<bool>("beamSpotCorr")),
  srcBeamSpot_(iConfig.getParameter<edm::InputTag>("srcBeamSpot")),
  genParticleSrc_(iConfig.getParameter<edm::InputTag>("genParticles")),
  qualityCutsPSet_(iConfig.getParameter<edm::ParameterSet>("qualityCuts"))
{

  edm::ParameterSet isolationQCuts = qualityCutsPSet_.getParameterSet("isolationQualityCuts");
  qcuts_.reset(new reco::tau::RecoTauQualityCuts(isolationQCuts));
  vertexAssociator_.reset(new reco::tau::RecoTauVertexAssociator(qualityCutsPSet_));
  std::pair<edm::ParameterSet, edm::ParameterSet> puFactorizedIsoQCuts = reco::tau::factorizePUQCuts(isolationQCuts);
  if ( iConfig.exists("deltaBetaPUTrackPtCutOverride") ) {
    puFactorizedIsoQCuts.second.addParameter<double>("minTrackPt", iConfig.getParameter<double>("deltaBetaPUTrackPtCutOverride"));
  } else {
    puFactorizedIsoQCuts.second.addParameter<double>("minTrackPt", isolationQCuts.getParameter<double>("minGammaEt"));
  }
  pileupQcutsPUTrackSelection_.reset(new reco::tau::RecoTauQualityCuts(puFactorizedIsoQCuts.first));
  pileupQcutsGeneralQCuts_.reset(new reco::tau::RecoTauQualityCuts(puFactorizedIsoQCuts.second));

  pfCandSrc_ = iConfig.getParameter<edm::InputTag>("particleFlowSrc");

}

TauBlock::~TauBlock() { }

void TauBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneTau = new TClonesArray("vhtm::Tau");
  tree->Branch("Tau", &cloneTau, 32000, 2);
  tree->Branch("nTau", &fnTau, "fnTau/I");
}

void TauBlock::beginEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // NB: The use of the PV in this context is necessitated by its use in
  // applying quality cuts to the different objects in the isolation cone
  // The vertex associator contains the logic to select the appropriate vertex
  // We need to pass it the event so it can load the vertices.
  vertexAssociator_->setEvent(iEvent);

}

void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTau->Clear();
  fnTau = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_inputTag, taus);

  this->beginEvent(iEvent, iSetup);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleSrc_, genParticles);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(srcBeamSpot_, beamSpotHandle);

  edm::Handle<reco::PFCandidateCollection> pfCandHandle_;
  iEvent.getByLabel(pfCandSrc_, pfCandHandle_);

  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  
  if (taus.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(_vtxInputTag, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    if (_beamSpotCorr) {
      iEvent.getByLabel(_beamSpotInputTag, beamSpot);
    }

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
                                               it != taus->end(); ++it) {
      if (fnTau == kMaxTau) {
	edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << fnTau;
	break;
      }
      tauB = new ((*cloneTau)[fnTau++]) vhtm::Tau();

      //New variables
      vdouble newVars = produceNewVars(it, muons, primaryVertices, beamSpotHandle, trackBuilderHandle);
      tauB->isOneProng = newVars[0];
      tauB->isThreeProng = newVars[1];
      tauB->IPxyLead = newVars[2];
      tauB->sigma_IPxyLead = newVars[3];
      tauB->IPxySigLead = newVars[4];
      tauB->IPzLead = newVars[5];
      tauB->distSV_PVRef = newVars[6];
      tauB->distSV_PVSig = newVars[13];

      AlgebraicVector3 u1 = AlgebraicVector3(1., 0., 0.);
      AlgebraicVector3 u2 = AlgebraicVector3(0., 1., 0.);
      AlgebraicVector3 u3 = AlgebraicVector3(0., 0., 1.);
      //float tauPx = it->px();
      //float tauPy = it->py();
      //float tauPz = it->pz();

      //AlgebraicVector3 pvRecoVertex;
      //AlgebraicVector3 svRecoVertex;
      //AlgebraicVector3 pvGenVertex;
      //AlgebraicVector3 svGenVertex;

      std::vector<double> vtxPos; 
      if(iEvent.id().run() == 1){

		vtxPos = GenSecondaryVertex(it, genParticles);
   		std::cout<<"GenPrimVtx "<<vtxPos[0]<<" "<<vtxPos[1]<<" "<<vtxPos[2]<<std::endl;
   		std::cout<<"GenSecVtx "<<vtxPos[3]<<" "<<vtxPos[4]<<" "<<vtxPos[5]<<std::endl;
      		std::cout<<"RecoPrimVtx "<<newVars[7]<<" "<<newVars[8]<<" "<<newVars[9]<<std::endl;
      		std::cout<<"RecoSecVtx "<<newVars[10]<<" "<<newVars[11]<<" "<<newVars[12]<<std::endl;

      		AlgebraicVector3 taudir = AlgebraicVector3(vtxPos[6], vtxPos[7], vtxPos[8]);

      		SVfit_namespace::compLocalCoordinates(taudir, u1, u2, u3);

		//pvGenVertex = AlgebraicVector3(newVars[0], newVars[1], newVars[2]);
		//pvRecoVertex = AlgebraicVector3(newVars[7], newVars[8], newVars[9]);

		//AlgebraicVector3 pvRecoVertexLocal = SVfit_namespace::transformToLocalCoordinatesPos(pvGenVertex, u1, u2, u3);
		//AlgebraicVector3 pvGenVertexLocal = SVfit_namespace::transformToLocalCoordinatesPos(pvRecoVertex, u1, u2, u3);

		tauB->diffPVx = newVars[7] - vtxPos[0];
		tauB->diffPVy = newVars[8] - vtxPos[1];
		tauB->diffPVz = newVars[9] - vtxPos[2];

		float resPVx = newVars[7] - vtxPos[0];
		float resPVy = newVars[8] - vtxPos[1];
		float resPVz = newVars[9] - vtxPos[2];
		AlgebraicVector3 resPV = AlgebraicVector3(resPVx, resPVy, resPVz);

		AlgebraicVector3 resPVLocal = SVfit_namespace::transformToLocalCoordinatesPos(resPV, u1, u2, u3);

		tauB->diffPVxLocal = resPVLocal[0];
		tauB->diffPVyLocal = resPVLocal[1];
		tauB->diffPVzLocal = resPVLocal[2];

		//std::cout<<newVars[9]<<" "<<vtxPos[2]<<" "<<(newVars[9] - vtxPos[2])<<std::endl;

		if(newVars[1]){

			//svGenVertex = AlgebraicVector3(newVars[3], newVars[4], newVars[5]);
			//svRecoVertex = AlgebraicVector3(newVars[10], newVars[11], newVars[12]);

			//AlgebraicVector3 svRecoVertexLocal = SVfit_namespace::transformToLocalCoordinatesPos(svGenVertex, u1, u2, u3);
			//AlgebraicVector3 svGenVertexLocal = SVfit_namespace::transformToLocalCoordinatesPos(svRecoVertex, u1, u2, u3);

			tauB->diffSVx = newVars[10] - vtxPos[3];
			tauB->diffSVy = newVars[11] - vtxPos[4];
			tauB->diffSVz = newVars[12] - vtxPos[5];

			float resSVx = newVars[10] - vtxPos[3];
			float resSVy = newVars[11] - vtxPos[4];
			float resSVz = newVars[12] - vtxPos[5];
			AlgebraicVector3 resSV = AlgebraicVector3(resSVx, resSVy, resSVz);

			AlgebraicVector3 resSVLocal = SVfit_namespace::transformToLocalCoordinatesPos(resSV, u1, u2, u3);

			tauB->diffSVxLocal = resSVLocal[0];
			tauB->diffSVyLocal = resSVLocal[1];
			tauB->diffSVzLocal = resSVLocal[2];
	
		}

      }

	reco::VertexRef pv = vertexAssociator_->associatedVertex(*(it->pfJetRef()));
	if (pv.isNonnull() && it->leadPFChargedHadrCand().isNonnull() && TMath::Abs(pv->z() - it->vertex().z()) < 0.2 && it->tauID("decayModeFinding") > 0.5){
      if (tb::isValidRef(it->leadPFChargedHadrCand()) && tb::isValidRef(it->leadPFChargedHadrCand()->trackRef())) {
	reco::TrackRef trk = it->leadPFChargedHadrCand()->trackRef();
		std::cout<<"ngul a mamet3 "<<pv->x()<<std::endl;
		qcuts_->setPV(pv);
		std::cout<<"ngul a mamet3"<<std::endl;
		qcuts_->setLeadTrack(it->leadPFChargedHadrCand());
		std::cout<<"ngul a mamet3"<<std::endl;
		pileupQcutsGeneralQCuts_->setPV(pv);
		pileupQcutsGeneralQCuts_->setLeadTrack(it->leadPFChargedHadrCand());
		pileupQcutsPUTrackSelection_->setPV(pv);
		pileupQcutsPUTrackSelection_->setLeadTrack(it->leadPFChargedHadrCand());
		std::cout<<"ngul a mamet3"<<std::endl;
		std::vector<double> isoQuantities = produceIso(it, pfCandHandle_, qcuts_, pileupQcutsPUTrackSelection_, pileupQcutsGeneralQCuts_, 0.5, 0.8, 0.3);}}

      // Store Tau variables
      tauB->eta    = it->eta();
      tauB->phi    = it->phi();
      tauB->pt     = it->pt();
      tauB->energy = it->energy();
      tauB->charge = it->charge();
      if (tb::isValidRef(it->leadPFChargedHadrCand()) && tb::isValidRef(it->leadPFChargedHadrCand()->trackRef())) {
	reco::TrackRef trk = it->leadPFChargedHadrCand()->trackRef();
        tauB->leadTrkPt      = trk->pt();
        tauB->leadTrkPtError = trk->ptError();
        tauB->leadTrkEta     = trk->eta();
        tauB->leadTrkPhi     = trk->phi();
        tauB->leadTrkCharge  = trk->charge();

        double trkd0 = trk->d0();
        double trkdz = trk->dz();
        if (_beamSpotCorr) {
          if (beamSpot.isValid()) {
            trkd0 = -(trk->dxy(beamSpot->position()));
            trkdz = trk->dz(beamSpot->position());
          }
          else
            edm::LogError("MuonsBlock") << "Error >> Failed to get BeamSpot for label: "
                                        << _beamSpotInputTag;
        }
        tauB->leadTrkD0      = trkd0;
        tauB->leadTrkD0Error = trk->d0Error();
        tauB->leadTrkDz      = trkdz;
        tauB->leadTrkDzError = trk->dzError();
        if (0) std::cout << trk->pt() << " "
                         << trk->eta() << " "
                         << trk->phi() << " "
                         << trk->charge() << " "
                         << trk->d0() <<  " "
                         << trk->dz() 
                         << std::endl;

        if (primaryVertices.isValid()) {
          edm::LogInfo("TauBlock") << "Total # Primary Vertices: " << primaryVertices->size();

          // IP of leadPFChargedHadrCand wrt event PV
          reco::VertexCollection::const_iterator vit = primaryVertices->begin(); // Highest sumPt vertex
          tauB->dxyPV = trk->dxy(vit->position());
          tauB->dzPV  = trk->dz(vit->position());

          // IP of leadPFChargedHadrCand wrt closest PV
          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDz = 9999.;
          double vertexDxy = 9999.;
          for (reco::VertexCollection::const_iterator vit  = primaryVertices->begin();
                                                      vit != primaryVertices->end(); ++vit) {
            double dxy = trk->dxy(vit->position());
            double dz  = trk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDxy = dxy;
              vertexDz = dz;
            }
          }
          tauB->vtxIndex = indexVtx;
          tauB->vtxDxy   = vertexDxy;
          tauB->vtxDz    = vertexDz;
          if (0) std::cout << indexVtx << " " << vertexDxy << " " << vertexDz << std::endl;
	}
        else {
  	  edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: "
                                    << _vtxInputTag;
        }
      }
      // Leading particle pT
      tauB->leadChargedParticlePt = it->leadPFChargedHadrCand().isNonnull() 
                                      ? it->leadPFChargedHadrCand()->pt(): 0.;
      tauB->leadNeutralParticlePt = it->leadPFNeutralCand().isNonnull() 
                                      ? it->leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = it->leadPFCand().isNonnull() 
                                      ? it->leadPFCand()->et(): 0.;      
      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = it->signalPFCands().size();
      
      tauB->numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = it->isolationPFCands().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = it->isolationPFGammaCandsEtSum();

      // Signal Constituents
      // Charged hadrons
      int indx = 0;
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFChargedHadrCands().begin(); 
    	                                              iCand != it->signalPFChargedHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          tauB->sigChHadCandPt[indx]  = cand.pt();
          tauB->sigChHadCandEta[indx] = cand.eta();
          tauB->sigChHadCandPhi[indx] = cand.phi();
        }
      }  
      // Neutral hadrons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFNeutrHadrCands().begin(); 
    	                                              iCand != it->signalPFNeutrHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->sigNeHadCandPt[indx]  = cand.pt();
          tauB->sigNeHadCandEta[indx] = cand.eta();
          tauB->sigNeHadCandPhi[indx] = cand.phi();
        }
      }  
      // Photons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFGammaCands().begin(); 
    	                                              iCand != it->signalPFGammaCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->sigGammaCandPt[indx]  = cand.pt();
          tauB->sigGammaCandEta[indx] = cand.eta();
          tauB->sigGammaCandPhi[indx] = cand.phi();
        }
      }  
      // Isolation Constituents
      // Charged hadrons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFChargedHadrCands().begin(); 
    	                                              iCand != it->isolationPFChargedHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          tauB->isoChHadCandPt[indx]  = cand.pt();
          tauB->isoChHadCandEta[indx] = cand.eta();
          tauB->isoChHadCandPhi[indx] = cand.phi();
        }
      }  
      
      // Neutral hadrons
      double ptSum = 0;
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFNeutrHadrCands().begin(); 
	                                              iCand != it->isolationPFNeutrHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->isoNeHadCandPt[indx]  = cand.pt();
          tauB->isoNeHadCandEta[indx] = cand.eta();
          tauB->isoNeHadCandPhi[indx] = cand.phi();
        }
	ptSum += std::abs(cand.pt());
      }  
      tauB->ptSumPFNeutralHadronsIsoCone = ptSum;

      // Photons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFGammaCands().begin(); 
      	                                              iCand != it->isolationPFGammaCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->isoGammaCandPt[indx]  = cand.pt();
          tauB->isoGammaCandEta[indx] = cand.eta();
          tauB->isoGammaCandPhi[indx] = cand.phi();
        }
      }  

      // tau id. discriminators
      tauB->decayModeFinding = it->tauID("decayModeFinding");

      // discriminators against muons
      tauB->againstMuonLoose      = it->tauID("againstMuonLoose");
      tauB->againstMuonMedium     = it->tauID("againstMuonMedium");
      tauB->againstMuonTight      = it->tauID("againstMuonTight");

      tauB->againstMuonLoose2     = it->tauID("againstMuonLoose2");
      tauB->againstMuonMedium2    = it->tauID("againstMuonMedium2");
      tauB->againstMuonTight2     = it->tauID("againstMuonTight2");

      // discriminators against electrons
      tauB->againstElectronLoose  = it->tauID("againstElectronLoose");
      tauB->againstElectronMedium = it->tauID("againstElectronMedium");
      tauB->againstElectronTight  = it->tauID("againstElectronTight");
      tauB->againstElectronMVA    = it->tauID("againstElectronMVA");

      tauB->againstElectronVLooseMVA2  = it->tauID("againstElectronVLooseMVA2");
      tauB->againstElectronLooseMVA2   = it->tauID("againstElectronLooseMVA2");
      tauB->againstElectronMediumMVA2  = it->tauID("againstElectronMediumMVA2");
      tauB->againstElectronTightMVA2   = it->tauID("againstElectronTightMVA2");
      
      tauB->againstElectronLooseMVA3  = it->tauID("againstElectronLooseMVA3");
      tauB->againstElectronMediumMVA3  = it->tauID("againstElectronMediumMVA3");
      tauB->againstElectronTightMVA3  = it->tauID("againstElectronTightMVA3");
      tauB->againstElectronVTightMVA3  = it->tauID("againstElectronVTightMVA3");

      // Obsolete
      tauB->pfElectronMVA         = it->leadPFCand().isNonnull() 
                                  ? it->leadPFCand()->mva_e_pi() : 1.;

      tauB->byVLooseIsolation  = it->tauID("byVLooseIsolation");
      tauB->byLooseIsolation   = it->tauID("byLooseIsolation");
      tauB->byMediumIsolation  = it->tauID("byMediumIsolation");
      tauB->byTightIsolation   = it->tauID("byTightIsolation");

      // MVA based isolation
      tauB->byLooseIsolationMVA                    = it->tauID("byLooseIsolationMVA"); 
      tauB->byMediumIsolationMVA                   = it->tauID("byMediumIsolationMVA"); 
      tauB->byTightIsolationMVA                    = it->tauID("byTightIsolationMVA"); 

      tauB->byLooseIsolationMVA2                   = it->tauID("byLooseIsolationMVA2"); 
      tauB->byMediumIsolationMVA2                  = it->tauID("byMediumIsolationMVA2"); 
      tauB->byTightIsolationMVA2                   = it->tauID("byTightIsolationMVA2"); 
  
      // DB Corrected Isolation
      tauB->byVLooseCombinedIsolationDeltaBetaCorr = it->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      tauB->byLooseCombinedIsolationDeltaBetaCorr  = it->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      tauB->byMediumCombinedIsolationDeltaBetaCorr = it->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      tauB->byTightCombinedIsolationDeltaBetaCorr  = it->tauID("byTightCombinedIsolationDeltaBetaCorr");

      tauB->byVLooseIsolationDeltaBetaCorr         = it->tauID("byVLooseIsolationDeltaBetaCorr");
      tauB->byLooseIsolationDeltaBetaCorr          = it->tauID("byLooseIsolationDeltaBetaCorr");
      tauB->byMediumIsolationDeltaBetaCorr         = it->tauID("byMediumIsolationDeltaBetaCorr");
      tauB->byTightIsolationDeltaBetaCorr          = it->tauID("byTightIsolationDeltaBetaCorr");

      tauB->byLooseCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      tauB->byMediumCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tauB->byTightCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");

      // kinematic variables for PFJet associated to PFTau
      if (tb::isValidRef(it->pfJetRef())) {
	const reco::PFJetRef &jtr = it->pfJetRef();
        tauB->jetPt  = jtr->pt();
        tauB->jetEta = jtr->eta();
        tauB->jetPhi = jtr->phi();
        if (0) std::cout << jtr->pt() << " " << jtr->eta() << " " << jtr->phi() << std::endl;
      }

      // NEW quantities
      tauB->emFraction              = it->emFraction(); 
      tauB->maximumHCALPFClusterEt  = it->maximumHCALPFClusterEt();
      tauB->ecalStripSumEOverPLead  = it->ecalStripSumEOverPLead();
      tauB->bremsRecoveryEOverPLead = it->bremsRecoveryEOverPLead();
      tauB->hcalTotOverPLead        = it->hcalTotOverPLead();
      tauB->hcalMaxOverPLead        = it->hcalMaxOverPLead();
      tauB->hcal3x3OverPLead        = it->hcal3x3OverPLead();

      tauB->etaetaMoment = it->etaetaMoment();
      tauB->phiphiMoment = it->phiphiMoment();
      tauB->phiphiMoment = it->etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      tauB->vx = vertex.x();             
      tauB->vy = vertex.y();             
      tauB->vz = vertex.z();             

      tauB->zvertex = it->vz(); // distance from the primary vertex
      tauB->mass    = it->p4().M();
      tauB->ltsipt  = TMath::Abs(it->leadPFChargedHadrCandsignedSipt());
    }
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: " 
                              << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
