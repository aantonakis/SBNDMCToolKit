#ifndef RECO_TRACK_INFO_H
#define RECO_TRACK_INFO_H

// LarSoft Includes
#include "lardataobj/RecoBase/Shower.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"


// Standard Library Includes
#include <iostream>
#include <fstream>

// // ROOT Includes
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"


namespace analysis {


class RecoTrackInfo {

public:
        RecoTrackInfo() {}
        ~RecoTrackInfo() {}

        unsigned int fEventID;
  	unsigned int fRun;
	unsigned int fSubRun;
        // Attributes of Reconstructed Showers
	Short_t  ntracks;                  // Number of reconstructed tracks
        std::vector<Int_t>  ftrkID;        // track ID
	std::vector<Int_t>  ftrkPFPID; // ID of the corresponding PFP
	std::vector<Int_t>  ftrkIsPrimary; // Is it Primary according to Reco
	std::vector<Int_t>  ftrkNDaughters; // How many daughter PFPs
	std::vector<Int_t>  ftrkParentPFPID; // The parent PFP ID
  	std::vector<float> ftrkLengths;
        std::vector<float> ftrkstartx;
        std::vector<float> ftrkstarty;
        std::vector<float> ftrkstartz;
        //std::vector<float> ftrkstartd;
        std::vector<float> ftrkendx;
        std::vector<float> ftrkendy;
        std::vector<float> ftrkendz;
	std::vector<float> ftrkmom;
        //std::vector<float> ftrkendd;
        std::vector<float> ftrktheta;
        std::vector<float> ftrkphi;
        std::vector<float> ftrkstartdcosx;
        std::vector<float> ftrkstartdcosy;
        std::vector<float> ftrkstartdcosz;
  	std::vector<std::vector<float>> ftrkdEdx;
  	std::vector<std::vector<float>> ftrkResidualRange;
	std::vector<float> ftrkcosmicscore;

        typedef std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle> > trkPfpMap;
	typedef std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>>::const_iterator trkPfpMapIt;
	trkPfpMap trackPFParticleMap;

        // call this at the beginning of each LArSoft Event
	void clearRecoTrackData() {
          ftrkID.clear();
	  ftrkPFPID.clear();
	  ftrkIsPrimary.clear();
	  ftrkNDaughters.clear();
	  ftrkParentPFPID.clear();
	  ftrkLengths.clear();
	  ftrkstartx.clear();
	  ftrkstarty.clear();
	  ftrkstartz.clear();
	  //ftrkstartd.clear();
	  ftrkendx.clear();
	  ftrkendy.clear();
	  ftrkendz.clear();
	  ftrkmom.clear();
	  //ftrkendd.clear();
	  ftrktheta.clear();
	  ftrkphi.clear();
	  ftrkstartdcosx.clear();
	  ftrkstartdcosy.clear();
	  ftrkstartdcosz.clear();
	  ftrkdEdx.clear();
	  ftrkResidualRange.clear();
	  ftrkcosmicscore.clear();

	  trackPFParticleMap.clear();

        } // end of clearRecoShowerData

        void setRecoTrackAddresses(TTree *tree,  std::string tracklabel) {
  	  std::cout << "tracklabel " << tracklabel << std::endl;
	  tree->Branch("EventID", &fEventID);
	  tree->Branch("Run", &fRun);	  
	  tree->Branch("SubRun", &fSubRun);
	  tree->Branch("trkID", &ftrkID);
	  tree->Branch("trkPFPID", &ftrkPFPID);
	  tree->Branch("trkIsPrimary", &ftrkIsPrimary);
	  tree->Branch("trkNDaughters", &ftrkNDaughters);
	  tree->Branch("trkParentPFPID", &ftrkParentPFPID);
	  tree->Branch("trkLengths", &ftrkLengths);        
          tree->Branch("trkstartx", &ftrkstartx);	  
          tree->Branch("trkstarty", &ftrkstarty);	  
          tree->Branch("trkstartz", &ftrkstartz);	  
	  //tree->Branch("trkstartd", &ftrkstartd);   
	  tree->Branch("trkendx", &ftrkendx);
	  tree->Branch("trkendy", &ftrkendy);
	  tree->Branch("trkendz", &ftrkendz);
	  tree->Branch("trkmom", &ftrkmom);
	  //tree->Branch("trkendd", &ftrkendd);
	  tree->Branch("trktheta", &ftrktheta);
	  tree->Branch("trkphi", &ftrkphi);
	  tree->Branch("trkstartdcosx", &ftrkstartdcosx);
	  tree->Branch("trkstartdcosy", &ftrkstartdcosy);
	  tree->Branch("trkstartdcosz", &ftrkstartdcosz);
	  tree->Branch("trkdEdx", &ftrkdEdx);
	  tree->Branch("trkResidualRange", &ftrkResidualRange);
	  tree->Branch("trkcosmicscore", &ftrkcosmicscore);

        } // end of setRecoShowerAddresses


        void makeTrkPFPMap(art::Event const& e, std::string tracklabel, art::Handle<std::vector<recob::PFParticle>> pfpHandle) {
      
	  art::FindManyP<recob::Track> ftrk(pfpHandle, e, tracklabel);
          
	  for(unsigned int i = 0; i < pfpHandle->size(); ++i) {
            art::Ptr< recob::PFParticle > pfp( pfpHandle, i );

            // Get Track association
            std::vector< art::Ptr<recob::Track> > trkAssn = ftrk.at(pfp->Self());
                 
            if(trkAssn.size()  > 1){
              std::cout << "PFParticle has " << trkAssn.size() << " associated tracks, should only have 1 or 0 ";
              continue;
            }
            if(trkAssn.size() == 0) continue;
            trackPFParticleMap.emplace(trkAssn[0],pfp);

          }
         
        } // end of makeShwPFPMap

        void fillRecoTrackTree(TTree *tree, art::Event const& e, std::string tracklabel, std::string pfplabel, std::string calolabel, std::string cosmictaggerlabel) {
        
	  // Want to record what event we are in
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();

	  // Load tracks from Pandora or a different producer
  	  art::Handle<std::vector<recob::Track>> trackHandle;
          std::vector<art::Ptr<recob::Track>> trackVec;
  	  if (e.getByLabel(tracklabel, trackHandle))
            art::fill_ptr_vector(trackVec, trackHandle);
          
	  if (trackVec.empty())
            return;

  	  // Get Calorimetry info for each track  
  	  art::FindManyP<anab::Calorimetry> trackCaloAssns(trackHandle, e, calolabel);
          
	  // Loop over all of the reco tracks in this event
    	  for (const art::Ptr<recob::Track> &trk : trackVec) {
  	    trkPfpMapIt it;
            it = trackPFParticleMap.find(trk);
            if(it == trackPFParticleMap.end()) continue; // check if it's in the trackPFP map
            art::Ptr<recob::PFParticle> pfpAssn = it->second;

	    // First, let's match the pfp to this track
	    ftrkPFPID.push_back(pfpAssn->Self());
	    ftrkIsPrimary.push_back(pfpAssn->IsPrimary());
	    ftrkNDaughters.push_back(pfpAssn->NumDaughters());
	    ftrkParentPFPID.push_back(pfpAssn->Parent());
            ftrkID.push_back(trk->ID()); 	    
            ftrkLengths.push_back(trk->Length());
    	
 	    // Not sure if this will work
	    // track momentum algorithm based on track range
	    trkf::TrackMomentumCalculator trkm;

	    int ntraj = trk->NumberTrajectoryPoints();
            if (ntraj > 0) {
              const auto& pos       = trk->Vertex();
              const auto& dir_start = trk->VertexDirection();
              //const auto& dir_end   = trk->EndDirection();
              const auto& end       = trk->End();
	      double mom = 0.;
              //tlen        = length(track);
              if(trk->HasMomentum() > 0)
                mom = trk->VertexMomentum();
          
              //double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
              //double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
              //double dpos = bdist(pos);
              //double dend = bdist(end);
        
              ftrkstartx.push_back(pos.X());
              ftrkstarty.push_back(pos.Y());
              ftrkstartz.push_back(pos.Z());
              //ftrkstartd.push_back(dpos);
              ftrkendx.push_back(end.X());
              ftrkendy.push_back(end.Y());
              ftrkendz.push_back(end.Z());
	      ftrkmom.push_back(mom);
              //ftrkendd.push_back(dend);
              ftrktheta.push_back(dir_start.Theta());
              ftrkphi.push_back(dir_start.Phi());
              ftrkstartdcosx.push_back(dir_start.X());
              ftrkstartdcosy.push_back(dir_start.Y());
              ftrkstartdcosz.push_back(dir_start.Z());

	    } // number of traj points

            // Calo Block

            const std::vector<art::Ptr<anab::Calorimetry>> trackCalos(trackCaloAssns.at(trk.key()));
            for (const art::Ptr<anab::Calorimetry> &calo: trackCalos) {
              const int planeNum(calo->PlaneID().Plane);
              // if it is not on the collection plane --> skip for now
              if (planeNum != 2)
                continue;
              ftrkdEdx.push_back(calo->dEdx());
	      ftrkResidualRange.push_back(calo->ResidualRange());
            }
                 
            // End Calo Block      

	    // Let's see if we can flag stuff as being a potential cosmic
	    art::FindManyP<anab::CosmicTag> fmct(trackHandle, e, cosmictaggerlabel);
	    if (fmct.isValid()) {
              ftrkcosmicscore.push_back(fmct.at(trk.key()).at(0)->CosmicScore());	    
	    }
	    else {
	      ftrkcosmicscore.push_back(-1.);
	      //std::cout << "Invalid Cosmic Tagger Association" << std::endl;
	    }
          } // end of loop over tracks

	  // TODO Note: we could put the Fill statement in the analyzer module if we decide to
	  // incorporate more info outside of this simple loop

          tree->Fill();
	 
        } // end of fillRecoTrackTree


private:
  // might not need this


};

} // end namespace analysis

#endif






