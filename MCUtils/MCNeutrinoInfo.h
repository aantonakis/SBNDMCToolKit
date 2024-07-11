#ifndef MC_NEUTRINO_INFO_H
#define MC_NEUTRINO_INFO_H

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
#include "nusimdata/SimulationBase/GTruth.h"

// Standard Library Includes
#include <iostream>
#include <fstream>

// // ROOT Includes
#include "TObject.h"
#include "TVector3.h"

// Custom includes
#include "PFParticleInfo.h"
#include "TruthMatcher.h"
#include "GeoHelper.h"
#include "MCParticleInfo.h"

namespace analysis {


class MCNeutrinoInfo {

public:
        MCNeutrinoInfo() {}
        ~MCNeutrinoInfo() {}


        unsigned int fEventID;
	unsigned int fRun;
	unsigned int fSubRun;	

	std::vector<Int_t> nu_ids;

	//Int_t                  mcevts_truth;    //number of neutrino Int_teractions in the spill
        std::vector<Int_t>     nuScatterCode_truth; //Scattering code given by Genie for each neutrino
        std::vector<Int_t>     nuID_truth;          //Unique ID of each true neutrino
        std::vector<Int_t>     nuPDG_truth;         //neutrino PDG code
        std::vector<Int_t>     ccnc_truth;          //0=CC 1=NC
        std::vector<Int_t>     mode_truth;          //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
        std::vector<Float_t>   enu_truth;           //true neutrino energy
        std::vector<Float_t>   Q2_truth;            //Momentum transfer squared
        std::vector<Float_t>   W_truth;             //hadronic invariant mass
        std::vector<Int_t>     hitnuc_truth;        //hit nucleon
        std::vector<Float_t>   nuvtxx_truth;        //neutrino vertex x
        std::vector<Float_t>   nuvtxy_truth;        //neutrino vertex y
        std::vector<Float_t>   nuvtxz_truth;        //neutrino vertex z
        std::vector<Float_t>   nu_dcosx_truth;      //neutrino dcos x
        std::vector<Float_t>   nu_dcosy_truth;      //neutrino dcos y
        std::vector<Float_t>   nu_dcosz_truth;      //neutrino dcos z
        std::vector<Float_t>   lep_mom_truth;       //lepton momentum
        std::vector<Float_t>   lep_dcosx_truth;     //lepton dcos x
        std::vector<Float_t>   lep_dcosy_truth;     //lepton dcos y
        std::vector<Float_t>   lep_dcosz_truth;     //lepton dcos z
	std::vector<Int_t>     best_slice_id;
	std::vector<float>     best_slice_purity;
	std::vector<float>     best_slice_completeness;
	std::vector<bool>      is_active_volume;   // is the neutrino interaction vertex in the active volume
        //flux information
        //std::vector<Float_t>   tpx_flux;        //Px of parent particle leaving BNB target
        //std::vector<Float_t>   tpy_flux;        //Py of parent particle leaving BNB target
        //std::vector<Float_t>   tpz_flux;        //Pz of parent particle leaving BNB target
        //std::vector<Int_t>     tptype_flux;     //Type of parent particle leaving BNB target

	// call this at the beginning of each LArSoft Event
	void clearMCNeutrinoData() {
	
	  nu_ids.clear();
          nuScatterCode_truth.clear(); //Scattering code given by Genie for each neutrino
          //nuID_truth.clear();          //Unique ID of each true neutrino
          nuPDG_truth.clear();         //neutrino PDG code
          ccnc_truth.clear();          //0=CC 1=NC
          mode_truth.clear();          //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
          enu_truth.clear();           //true neutrino energy
          Q2_truth.clear();            //Momentum transfer squared
          W_truth.clear();             //hadronic invariant mass
          hitnuc_truth.clear();        //hit nucleon
          nuvtxx_truth.clear();        //neutrino vertex x
          nuvtxy_truth.clear();        //neutrino vertex y
          nuvtxz_truth.clear();        //neutrino vertex z
          nu_dcosx_truth.clear();      //neutrino dcos x
          nu_dcosy_truth.clear();      //neutrino dcos y
          nu_dcosz_truth.clear();      //neutrino dcos z
          lep_mom_truth.clear();       //lepton momentum
          lep_dcosx_truth.clear();     //lepton dcos x
          lep_dcosy_truth.clear();     //lepton dcos y
          lep_dcosz_truth.clear();     //lepton dcos z
	
	  best_slice_id.clear();
	  best_slice_purity.clear();
	  best_slice_completeness.clear();

	  is_active_volume.clear(); 

	}

        void setMCNeutrinoAddresses(TTree *tree) {

	  tree->Branch("EventID", &fEventID);
	  tree->Branch("Run", &fRun);
	  tree->Branch("SubRun", &fSubRun);
	  tree->Branch("nu_ids", &nu_ids);
	  tree->Branch("nuScatterCode_truth", &nuScatterCode_truth);
	  tree->Branch("nuPDG_truth", &nuPDG_truth);
	  tree->Branch("ccnc_truth", &ccnc_truth);
	  tree->Branch("mode_truth", &mode_truth);
	  tree->Branch("enu_truth", &enu_truth);
	  tree->Branch("Q2_truth", &Q2_truth);
	  tree->Branch("W_truth", &W_truth);
	  tree->Branch("hitnuc_truth", &hitnuc_truth);
	  tree->Branch("nuvtxx_truth", &nuvtxx_truth);
	  tree->Branch("nuvtxy_truth", &nuvtxy_truth);
	  tree->Branch("nuvtxz_truth", &nuvtxz_truth);
	  tree->Branch("nu_dcosx_truth", &nu_dcosx_truth);
	  tree->Branch("nu_dcosy_truth", &nu_dcosy_truth);
	  tree->Branch("nu_dcosz_truth", &nu_dcosz_truth);
	  tree->Branch("lep_mom_truth", &lep_mom_truth);
	  tree->Branch("lep_dcosx_truth", &lep_dcosx_truth);
	  tree->Branch("lep_dcosy_truth", &lep_dcosx_truth);
	  tree->Branch("lep_dcosz_truth", &lep_dcosx_truth);
	  tree->Branch("best_slice_id", &best_slice_id);
	  tree->Branch("best_slice_purity", &best_slice_purity);
	  tree->Branch("best_slice_completeness", &best_slice_completeness);
	  tree->Branch("is_active_volume", &is_active_volume);

	}

	// Make sure to call this after the PFParticleInfo fill function has been called --> we can use that info directly instead of using art stuff

	// Find the Reco slice that best matches a single true neutrino
	void findSlice(art::Event const& e, art::Ptr<simb::MCTruth>& mctruth, std::string hitlabel, std::string pfplabel, detinfo::DetectorClocksData const& clockData) {
	  
	  // Custom object to help us with truth matching
	  TruthMatcher truthMatch;
	 
	   // Get the hits in this event and check which ones belong to the neutrino
 	  art::Handle< std::vector<recob::Hit> > hitListHandle;
  	  std::vector<art::Ptr<recob::Hit> > hitlist;
  	  if (e.getByLabel(hitlabel, hitListHandle))
    	    art::fill_ptr_vector(hitlist, hitListHandle);

 	  // Get total number of hits for this neutrino
	  int nhits_originNu_total = 0;
          for(auto const hit : hitlist){
            Int_t trkId;
            if( truthMatch.HitTruthId(clockData, hit,trkId) ) {
              //art::Ptr<simb::MCTruth> mctruth;
              if( truthMatch.TrackMatchIdToMCTruth(trkId, mctruth) ) {
                nhits_originNu_total++;
              }
            }
          }

	  // Get Slices in this event
    	  art::Handle< std::vector<recob::Slice> > sliceHandle;
    	  std::vector< art::Ptr<recob::Slice> > slicelist;
    	  if (e.getByLabel(pfplabel,sliceHandle))
      	    art::fill_ptr_vector(slicelist, sliceHandle);
	  
	  // Get Assns bw hits and slices
	  art::FindManyP<recob::Hit> fm_slicehits(sliceHandle, e, pfplabel);
	  
	  float best_nu_score = 0.; // For now, let's make this purity + completeness 
	  float best_purity = 0.;
	  float best_completeness = 0.;
	  Int_t best_id = -1; // TODO check that this is indeed a nonsensical id --> there might actually be an id of -1

	  // Loop over slices and find the one that matches the neutrino the best
  	  for (const art::Ptr<recob::Slice> &slc: slicelist) {
	    std::vector<art::Ptr<recob::Hit>> sliceHits = fm_slicehits.at(slc.key());

	    int nhits_slice = (int)sliceHits.size();
	    // Initialize the current purity and completeness for this slice         
	    float purity = 0.;
	    float completeness = 0.;
	    int nhits_originNu = 0;

	    // Loop over the hits in this slice
	    for(auto const hit : sliceHits){
            
	      Int_t trkId;
              if( truthMatch.HitTruthId(clockData,hit,trkId) ) {
                if( truthMatch.TrackMatchIdToMCTruth(trkId, mctruth) ) {
                    nhits_originNu++;
	        }
	      }
	    } // End loop over hits

	    purity += float(nhits_originNu);
	    completeness += float(nhits_originNu);
	    
	    if (purity < 1.)
	      continue;
	    
	    purity /= float(nhits_slice);
	    completeness /= float(nhits_originNu_total);
	
	    float temp_score = purity * completeness;
	    if (temp_score > best_nu_score) {
	      best_nu_score = temp_score;
	      best_purity = purity;
	      best_completeness = completeness;
	      best_id = slc->ID();
	    }	

	  } // End loop over slices in Event

	  best_slice_id.push_back(best_id);
	  best_slice_purity.push_back(best_purity);
	  best_slice_completeness.push_back(best_completeness);

	} // end of findSlice


        void fillMCNeutrinoTree(TTree *tree, art::Event const& e, std::string nugenlabel, std::string hitlabel, 
					                          std::string pfplabel, detinfo::DetectorClocksData const& clockData) {

	  
  	  //MCParticleInfo mcNuPartInfo;

	  std::cout << "In custom MCNeutrino Fill function" << std::endl;        
	  // Want to record what event we are in
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();


	  // Custom object to help us with truth matching
	  TruthMatcher truthMatch;
	  
	  // * MC truth information --> How we will get Neutrino Info
  	  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  	  std::vector<art::Ptr<simb::MCTruth> > mclist;
  	  if (e.getByLabel(nugenlabel, mctruthListHandle))
  	    art::fill_ptr_vector(mclist, mctruthListHandle);
	
	  // Get GTruth Info for a given scattering code
	  art::FindManyP< simb::GTruth > fmgt( mctruthListHandle, e, nugenlabel );

	  unsigned int nu_id = 0;

	  // Loop over the MC truths and find the neutrinos
	  for (unsigned int i_mctruth = 0; i_mctruth < mclist.size(); i_mctruth++){
            // fetch a mctruth
            art::Ptr<simb::MCTruth> curr_mctruth = mclist[i_mctruth];
            //Check if it's a neutrino
            if (!curr_mctruth->NeutrinoSet()) continue;
	    nu_ids.push_back(nu_id);

	    // Genie Truth associated only with the neutrino
	    std::vector< art::Ptr<simb::GTruth> > mcgtAssn = fmgt.at(i_mctruth);

	    nuScatterCode_truth.push_back(mcgtAssn[0]->fGscatter);

            nuPDG_truth.push_back(curr_mctruth->GetNeutrino().Nu().PdgCode());
            ccnc_truth.push_back(curr_mctruth->GetNeutrino().CCNC());
            mode_truth.push_back(curr_mctruth->GetNeutrino().Mode());
            Q2_truth.push_back(curr_mctruth->GetNeutrino().QSqr());
            W_truth.push_back(curr_mctruth->GetNeutrino().W());
            hitnuc_truth.push_back(curr_mctruth->GetNeutrino().HitNuc());
            enu_truth.push_back(curr_mctruth->GetNeutrino().Nu().E());
            nuvtxx_truth.push_back(curr_mctruth->GetNeutrino().Nu().Vx());
            nuvtxy_truth.push_back(curr_mctruth->GetNeutrino().Nu().Vy());
            nuvtxz_truth.push_back(curr_mctruth->GetNeutrino().Nu().Vz());

	    // Check if the vertex is in the active volume
	    GeoHelper geo_help;
	    if ( geo_help.isActive(curr_mctruth->GetNeutrino().Nu().Vx(), curr_mctruth->GetNeutrino().Nu().Vy(), curr_mctruth->GetNeutrino().Nu().Vz()) ) {
	      is_active_volume.push_back(true);
	    }
	    else {
	      is_active_volume.push_back(false);
	    }
            if (curr_mctruth->GetNeutrino().Nu().P()){
              nu_dcosx_truth.push_back(curr_mctruth->GetNeutrino().Nu().Px()/curr_mctruth->GetNeutrino().Nu().P());
              nu_dcosy_truth.push_back(curr_mctruth->GetNeutrino().Nu().Py()/curr_mctruth->GetNeutrino().Nu().P());
              nu_dcosz_truth.push_back(curr_mctruth->GetNeutrino().Nu().Pz()/curr_mctruth->GetNeutrino().Nu().P());
            }
	    else {
	      // Put garbage values
	      nu_dcosx_truth.push_back(-9999.);
	      nu_dcosy_truth.push_back(-9999.);
	      nu_dcosz_truth.push_back(-9999.);
	    }
            lep_mom_truth.push_back(curr_mctruth->GetNeutrino().Lepton().P());
            if (curr_mctruth->GetNeutrino().Lepton().P()){
              lep_dcosx_truth.push_back(curr_mctruth->GetNeutrino().Lepton().Px()/curr_mctruth->GetNeutrino().Lepton().P());
              lep_dcosy_truth.push_back(curr_mctruth->GetNeutrino().Lepton().Py()/curr_mctruth->GetNeutrino().Lepton().P());
              lep_dcosz_truth.push_back(curr_mctruth->GetNeutrino().Lepton().Pz()/curr_mctruth->GetNeutrino().Lepton().P());
            }

	    // Finished Filling truth info for this neutrino
	    // Now, Let's try to find the best slice for this neutrino
	    findSlice(e, curr_mctruth, hitlabel, pfplabel, clockData);
	   
	    nu_id += 1;

	  }// end loop over mctruths

	  tree->Fill();

	} // end of fill function


private:


};

} // end namespace analysis

#endif


