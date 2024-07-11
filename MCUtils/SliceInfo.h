#ifndef SLICE_INFO_H
#define SLICE_INFO_H

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


class SliceInfo {

public:
        SliceInfo() {}
        ~SliceInfo() {}

        unsigned int fEventID;
  	unsigned int fRun;
	unsigned int fSubRun;
        // Attributes of PFParticles
        
	std::vector<std::vector<Int_t>>   fpfpID;        // pfp ID
        std::vector<std::vector<Int_t>>   fpfpParentID;  // pfp Parent ID
        std::vector<Int_t>  fpfpSliceID;   // Slice ID for this pfp
	std::vector<Float_t> fOpT0Slice;
	std::vector<std::vector<bool>>    fpfpIsPrimary;
	std::vector<std::vector<float>>   fpfpNuScore;   // How neutrino-like is this pfp
	std::vector<std::vector<int>>     fpfpIsNu;      // Is it a Neutrino according to Pandora
	std::vector<std::vector<int>>     fpfpPDGCode;
        std::vector<std::vector<Float_t>> fpfpVtxX;
        std::vector<std::vector<Float_t>> fpfpVtxY;
        std::vector<std::vector<Float_t>> fpfpVtxZ;
	std::vector<std::vector<Float_t>> fpfpTrackScore;
	std::vector<std::vector<bool>>    fpfpIsClearCosmic;

	// call this at the beginning of each LArSoft Event
	void clearSliceData() {

          fpfpID.clear();
	  fpfpParentID.clear();
	  fpfpSliceID.clear();
	  fOpT0Slice.clear();
	  fpfpIsPrimary.clear();
	  fpfpNuScore.clear();
	  fpfpIsNu.clear();
	  fpfpPDGCode.clear();
	  fpfpVtxX.clear();
	  fpfpVtxY.clear();
	  fpfpVtxZ.clear();
          fpfpTrackScore.clear();
	  fpfpIsClearCosmic.clear();
	}

        void setSliceAddresses(TTree *tree,  std::string pfplabel) {
          std::string PFPLabel = pfplabel;
 	  std::string BranchName = "pfp_"+PFPLabel; 
  	  std::cout << "BranchName " << BranchName << std::endl;
	  tree->Branch("EventID", &fEventID);
	  tree->Branch("Run", &fRun);	  
	  tree->Branch("SubRun", &fSubRun);
	  tree->Branch("pfpSliceID", &fpfpSliceID);
	  tree->Branch("pfpID", &fpfpID);
	  tree->Branch("pfpParentID", &fpfpParentID);
	  tree->Branch("OpT0Slice", &fOpT0Slice);
	  tree->Branch("pfpIsPrimary", &fpfpIsPrimary);
	  tree->Branch("pfpNuScore", &fpfpNuScore);
	  tree->Branch("pfpPDGCode", &fpfpPDGCode);
	  tree->Branch("pfpVtxX", &fpfpVtxX);
	  tree->Branch("pfpVtxY", &fpfpVtxY);
	  tree->Branch("pfpVtxZ", &fpfpVtxZ);
	  tree->Branch("pfpTrackScore", &fpfpTrackScore);
          tree->Branch("pfpIsClearCosmic", &fpfpIsClearCosmic);
	}

        void fillSliceTree(TTree *tree, art::Event const& e, std::string pfplabel, art::Handle<std::vector<recob::PFParticle>> pfpHandle, 
								  std::string vertexlabel, std::string opt0finderlabel) {
	  //std::cout << "In custom PFPInfo Fill function" << std::endl;        
	  // Want to record what event we are in
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();

	  // Get the Slices from this event
 	  art::Handle< std::vector<recob::Slice> > sliceHandle;
          std::vector< art::Ptr<recob::Slice> > slicelist;
          if (e.getByLabel(pfplabel, sliceHandle))
            art::fill_ptr_vector(slicelist, sliceHandle);

          art::FindManyP<recob::PFParticle> fm_slicepfp(slicelist, e, pfplabel);
          // Map PFPs to their IDs for retrieval of nu scores
          art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(pfpHandle, e, pfplabel);

	  art::FindManyP<recob::Vertex> fvtxPFPAssns(pfpHandle, e, vertexlabel);
	  //std::cout << "opT0FinderLabel " << opt0finderlabel << std::endl; 
	  art::FindManyP<sbn::OpT0Finder> sliceOpT0Assoc(sliceHandle, e, opt0finderlabel);

	  // Loop over slices and neutrinos
          for(auto const &slice : slicelist ) {
	    
	    Float_t fOpT0Slice_temp = -1.;
	    if (sliceOpT0Assoc.isValid()) {
	      std::cout << "sliceOpT0Assoc was VALID YES!!!!" << std::endl;
	      // timing info for slice
	      std::vector<art::Ptr<sbn::OpT0Finder>> opT0s = sliceOpT0Assoc.at(slice.key());
	      // Get the time with the best score 
	      std::sort(opT0s.begin(), opT0s.end(), [](auto const& a, auto const& b) {return a->score > b->score;});
	      //Float_t fOpT0Slice_temp = -1.;
	      if (opT0s.size() != 0) {
                fOpT0Slice_temp = opT0s[0]->time;
	      }

	    }
	    // We will fill the opT0 result during the PFP loop
	    fpfpSliceID.push_back(slice->ID());
	    fOpT0Slice.push_back(fOpT0Slice_temp);

            std::vector<art::Ptr<recob::PFParticle>> slicePFPs = fm_slicepfp.at(slice.key());
            
	    std::vector<Int_t>   pfpID_temp;        // pfp ID
            std::vector<Int_t>   pfpParentID_temp;  // pfp Parent ID
	    std::vector<bool>    pfpIsPrimary_temp;
	    std::vector<float>   pfpNuScore_temp;   // How neutrino-like is this pfp
	    std::vector<int>     pfpIsNu_temp;      // Is it a Neutrino according to Pandora
	    std::vector<int>     pfpPDGCode_temp;
            std::vector<Float_t> pfpVtxX_temp;
            std::vector<Float_t> pfpVtxY_temp;
            std::vector<Float_t> pfpVtxZ_temp;
	    std::vector<Float_t> pfpTrackScore_temp;
	    std::vector<bool>    pfpIsClearCosmic_temp;

	    for (auto const &pfp : slicePFPs) {
	      //std::cout << "Starting a new pfp" << std::endl;
	      pfpID_temp.push_back(pfp->Self());
	      pfpIsPrimary_temp.push_back(pfp->IsPrimary());
	      pfpParentID_temp.push_back(pfp->Parent());
	      pfpPDGCode_temp.push_back(pfp->PdgCode());
              if ((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)) { // look for neutrino
	        pfpIsNu_temp.push_back(1);              
	      }//endif neutrino
	      else {
	        pfpIsNu_temp.push_back(0);              
	      }
 	      // Loop over the MetaData
 	      std::cout << "pfp self " << pfp->Self() << std::endl;
	      std::cout << "pfp key" << pfp.key() << std::endl;
	      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = fm_pfpmd.at(pfp->Self()); 
	      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVeck = fm_pfpmd.at(pfp.key()); // Let's compare
	      //std::cout << "Starting MetaData loop" <<std::endl;
	      //std::cout << "Length of pfpMetaVec " << pfpMetaVec.size() << std::endl; 
	      //std::cout << "Length of pfpMetaVeck " << pfpMetaVeck.size() << std::endl; 
	      
	      if (pfpMetaVec.size() > 1) {
	        std::cout << "More than one MetaData for this pfp ???" << std::endl;
		pfpNuScore_temp.push_back(-1.);
		pfpTrackScore_temp.push_back(-1.);
		pfpIsClearCosmic_temp.push_back(false);
	 	//continue;
	      }

              for (auto const pfpMeta : pfpMetaVec) {
		if (pfpMetaVec.size() > 1)
		  continue;

		//std::cout << "Getting the PropertiesMap" << std::endl;
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
		//std::cout << "Getting the NuScore" << std::endl;
		try {
                  pfpNuScore_temp.push_back(propertiesMap.at("NuScore"));
		}
		catch (const std::exception& err) {
		  //std::cerr << "Error: " << err.what() << std::endl;
                  pfpNuScore_temp.push_back(-1.); // This is most likely clear cosmic or pfp that is not the neutrino  TODO
		}
		try {
		  pfpTrackScore_temp.push_back(propertiesMap.at("TrackScore"));
		}
		catch (const std::exception& err) {
		  pfpTrackScore_temp.push_back(-1.);
		}
		try {
		  pfpIsClearCosmic_temp.push_back(propertiesMap.at("IsClearCosmic"));
		}
		catch (const std::exception& err) {
		  pfpIsClearCosmic_temp.push_back(false);
		}

              }// End loop over metadata

	      // CafMaker Code
	      //auto const &propertiesMap (pfpMeta->GetPropertiesMap());
      	      //auto const &pfpTrackScoreIter(propertiesMap.find("TrackScore"));
              //srpfp.trackScore = (pfpTrackScoreIter == propertiesMap.end()) ? -5.f : pfpTrackScoreIter->second;


	      //std::cout << "About to get vertex associated with this pfp" << std::endl;
	      // Get vertex associated with this PFP
	      std::vector<art::Ptr<recob::Vertex>> vtxAssn = fvtxPFPAssns.at(pfp.key());
	      if (vtxAssn.size() == 1) {
		// We have a vertex
		//std::cout << "About to fill the Vertex of a pfp" << std::endl;
	        //std::cout << "X " << vtxAssn.at(0)->position().X() << std::endl;
	        //std::cout << "Y " << vtxAssn.at(0)->position().Y() << std::endl;
	        //std::cout << "Z " << vtxAssn.at(0)->position().Z() << std::endl;
	 	pfpVtxX_temp.push_back(vtxAssn.at(0)->position().X());
	 	pfpVtxY_temp.push_back(vtxAssn.at(0)->position().Y());
	 	pfpVtxZ_temp.push_back(vtxAssn.at(0)->position().Z());
              }
	      else {
		std::cout << "No Valid vertex for this pfp --> Give dummy values" << std::endl;
	 	pfpVtxX_temp.push_back(9999.);
	 	pfpVtxY_temp.push_back(9999.);
	 	pfpVtxZ_temp.push_back(9999.);

	      }

            }//endloop over slicePFPs
 
	    fpfpID.push_back(pfpID_temp);        // pfp ID
            fpfpParentID.push_back(pfpParentID_temp);  // pfp Parent ID
	    fpfpIsPrimary.push_back(pfpIsPrimary_temp);
	    fpfpNuScore.push_back(pfpNuScore_temp);   // How neutrino-like is this pfp
	    fpfpIsNu.push_back(pfpIsNu_temp);      // Is it a Neutrino according to Pandora
	    fpfpPDGCode.push_back(pfpPDGCode_temp);
            fpfpVtxX.push_back(pfpVtxX_temp);
            fpfpVtxY.push_back(pfpVtxY_temp);
            fpfpVtxZ.push_back(pfpVtxZ_temp);
	    fpfpTrackScore.push_back(pfpTrackScore_temp);	
	    fpfpIsClearCosmic.push_back(pfpIsClearCosmic_temp);

	    //pfpID_temp.clear();
	    //pfpParentID_temp.clear();
	    //pfpNuScore_temp.clear();
	    //pfpIsNu_temp.clear();
	    //pfpPDGCode_temp.clear();
	    //pfpVtxX_temp.clear();
	    //pfpVtxY_temp.clear();	    
	    //pfpVtxZ_temp.clear();

          }//endloop over slices

	  tree->Fill();

	} // end of fillSliceTree



private:
// might not need this



};

} // end namespace analysis

#endif



