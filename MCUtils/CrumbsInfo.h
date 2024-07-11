#ifndef CRUMBS_INFO_H
#define CRUMBS_INFO_H

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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"


// Standard Library Includes
#include <iostream>
#include <fstream>

// // ROOT Includes
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"


namespace analysis {


class CrumbsInfo {

public:
        CrumbsInfo() {}
        ~CrumbsInfo() {}

        unsigned int fEventID;
  	unsigned int fRun;
	unsigned int fSubRun;
        // Attributes of PFParticles
        std::vector<Int_t>  fcrumbsSliceID;   // Slice ID for this pfp
	std::vector<Float_t> fcrumbsScore;
	std::vector<Float_t> fcrumbsCCnumuScore;
	std::vector<Float_t> fcrumbsCCnueScore;
	std::vector<Float_t> fcrumbsNCScore;
	std::vector<Float_t> fcrumbsBestScore;
	std::vector<Int_t> fcrumbsBestId;

	// call this at the beginning of each LArSoft Event
	void clearCrumbsData() {
	  
	  fcrumbsSliceID.clear();
	  fcrumbsScore.clear();
      	  fcrumbsCCnumuScore.clear();
      	  fcrumbsCCnueScore.clear();
      	  fcrumbsNCScore.clear();
      	  fcrumbsBestScore.clear();
     	  fcrumbsBestId.clear();
	}

        void setCrumbsAddresses(TTree *tree) {
	  tree->Branch("EventID", &fEventID);
	  tree->Branch("Run", fRun);	  
	  tree->Branch("SubRun", &fSubRun);
	  tree->Branch("crumbsSliceID", &fcrumbsSliceID);
	  tree->Branch("fcrumbsScore", &fcrumbsScore);
	  tree->Branch("fcrumbsCCnumuScore", &fcrumbsCCnumuScore);
          tree->Branch("fcrumbsCCnueScore", &fcrumbsCCnueScore);
	  tree->Branch("fcrumbsNCScore", &fcrumbsNCScore);
	  tree->Branch("fcrumbsBestScore", &fcrumbsBestScore);	
  	  tree->Branch("fcrumbsBestId", &fcrumbsBestId);
	}

	// Fill this info Slice by Slice

        void fillCrumbsTree(TTree *tree, art::Event const& e, std::string pfplabel, std::string crumbslabel) {
	  //std::cout << "In custom Crumbs Fill function" << std::endl;        
	  // Want to record what event we are in
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();
	  
	  // Get the Slices from this event
 	  art::Handle< std::vector<recob::Slice> > sliceHandle;
          std::vector< art::Ptr<recob::Slice> > slicelist;
          if (e.getByLabel(pfplabel, sliceHandle))
            art::fill_ptr_vector(slicelist, sliceHandle);

	  art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssoc(sliceHandle, e, crumbslabel);
	  // Loop over slices and neutrinos
          for(auto const &slice : slicelist ) {

            fcrumbsSliceID.push_back(slice->ID());
	    if (sliceCrumbsAssoc.isValid()) {
	      //std::cout << "Found some CRUMBS" << std::endl;
	      //std::cout << "Try to access crumbs for slice " << slice.key() << std::endl;
	      // Crumbs info for slice
	      std::vector<art::Ptr<sbn::CRUMBSResult>> crumbs = sliceCrumbsAssoc.at(slice.key());
	      //std::cout << "Got the crumbs --> check size" << std::endl;
	      if (crumbs.size() != 1) {
	        std::cout << "More than 1 Crumbs Result for this slice --> SliceId " << slice->ID() << std::endl;
	       	fcrumbsScore.push_back(-1.);
      	        fcrumbsCCnumuScore.push_back(-1.);
      	        fcrumbsCCnueScore.push_back(-1.);
      	        fcrumbsNCScore.push_back(-1.);
      	        fcrumbsBestScore.push_back(-1.);
     	        fcrumbsBestId.push_back(-1);
	        continue; 
	      }
	      fcrumbsScore.push_back(crumbs[0]->score);
      	      fcrumbsCCnumuScore.push_back(crumbs[0]->ccnumuscore);
      	      fcrumbsCCnueScore.push_back(crumbs[0]->ccnuescore);
      	      fcrumbsNCScore.push_back(crumbs[0]->ncscore);
      	      fcrumbsBestScore.push_back(crumbs[0]->bestscore);
     	      fcrumbsBestId.push_back(crumbs[0]->bestid);
      
	    } // check for the crumbs
	    else {
	      // didn't find crumbs ;(
	      fcrumbsScore.push_back(-1.);
      	      fcrumbsCCnumuScore.push_back(-1.);
      	      fcrumbsCCnueScore.push_back(-1.);
      	      fcrumbsNCScore.push_back(-1.);
      	      fcrumbsBestScore.push_back(-1.);
     	      fcrumbsBestId.push_back(-1);

	    }
	

/*
 *      Other Stuff
	slice.crumbs_result.tpc.crlongtrackhitfrac = crumbs->tpc_CRFracHitsInLongestTrack;
      	slice.crumbs_result.tpc.crlongtrackdefl = crumbs->tpc_CRLongestTrackDeflection;
      	slice.crumbs_result.tpc.crlongtrackdiry = crumbs->tpc_CRLongestTrackDirY;
      	slice.crumbs_result.tpc.crnhitsmax = crumbs->tpc_CRNHitsMax;
      	slice.crumbs_result.tpc.nusphereeigenratio = crumbs->tpc_NuEigenRatioInSphere;
      	slice.crumbs_result.tpc.nufinalstatepfos = crumbs->tpc_NuNFinalStatePfos;
      	slice.crumbs_result.tpc.nutotalhits = crumbs->tpc_NuNHitsTotal;
      	slice.crumbs_result.tpc.nuspherespacepoints = crumbs->tpc_NuNSpacePointsInSphere;
      	slice.crumbs_result.tpc.nuvertexy = crumbs->tpc_NuVertexY;
      	slice.crumbs_result.tpc.nuwgtdirz = crumbs->tpc_NuWeightedDirZ;
      	slice.crumbs_result.tpc.stoppingchi2ratio = crumbs->tpc_StoppingChi2CosmicRatio;
      	slice.crumbs_result.pds.fmtotalscore = crumbs->pds_FMTotalScore;
      	slice.crumbs_result.pds.fmpe = crumbs->pds_FMPE;
      	slice.crumbs_result.pds.fmtime = crumbs->pds_FMTime;
      	slice.crumbs_result.pds.opt0score = crumbs->pds_OpT0Score;
      	slice.crumbs_result.pds.opt0measuredpe = crumbs->pds_OpT0MeasuredPE;
      	slice.crumbs_result.crt.trackscore = crumbs->crt_TrackScore;
      	slice.crumbs_result.crt.spscore = crumbs->crt_SPScore;
      	slice.crumbs_result.crt.tracktime = crumbs->crt_TrackTime;
      	slice.crumbs_result.crt.sptime = crumbs->crt_SPTime;

*/


          }//endloop over slices


	  tree->Fill();

	} // end of fillCrumbsTree
	


private:
// might not need this

};

} // end namespace analysis

#endif



