#ifndef TRUTH_MATCHER_H
#define TRUTH_MATCHER_H

// Standard Library Includes
#include <iostream>
#include <fstream>

// // ROOT Includes
#include "TObject.h"
#include "TVector3.h"


namespace analysis {


class TruthMatcher {

public:
        TruthMatcher() {}
        ~TruthMatcher() {}
	
	// BackTrackerService
	art::ServiceHandle<cheat::BackTrackerService> bt_serv;
     
	// This function calculates the leading MCParticle ID contributing to a hit and the
	// fraction of that hit's energy comes from that particle.

	void HitTruth( detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, 
				    Int_t& truthid, Float_t& truthidEnergyFrac, Float_t& energy, Float_t& numElectrons){


	  std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
  	  if( !trackIDEs.size() ) return;
 	  // Loop through and find the leading TrackIDE, and keep
 	  // track of the total energy of ALL IDEs.
	  Float_t maxe = 0;
  	  Float_t bestfrac = 0;
  	  Int_t bestid = 0;
  	  Int_t ne = 0;
  	  for(size_t i = 0; i < trackIDEs.size(); ++i){
    	    ne += trackIDEs[i].numElectrons;
    	    if( trackIDEs[i].energy > maxe ) {
      	      maxe = trackIDEs[i].energy;
      	      bestfrac = trackIDEs[i].energyFrac;
      	      bestid = trackIDEs[i].trackID;
    	    }
  	  }
  	  truthid = bestid;
  	  truthidEnergyFrac = bestfrac;
  	  energy = maxe;
  	  numElectrons = ne;

	} // end of HitTruth


	// helper function to get TrackID (returns false if no match found)
	bool HitTruthId( detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& mcid) {
  	  mcid = std::numeric_limits<int>::lowest();
  	  Float_t dummy1;
  	  Float_t dummy2;
  	  Float_t dummy3;
  	  HitTruth(clockData,hit,mcid,dummy1,dummy2,dummy3);
  	  if( mcid > std::numeric_limits<int>::lowest() ) return true;
  	  else return false;

	} // End of HitTruthID


	// Get MCTruth associated with TrackID using a try bracket to avoid
	// fatal exceptions (return false if no match or exception thrown)
	bool TrackIdToMCTruth( Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth ){
	  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	  bool matchFound = false;
	  try {
	    mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
	    matchFound = true;
	  } catch(...) {
	    std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
	    std::cout<<"Exception inside TrackIdToMCTruth " << std::endl;;
	    matchFound = false;
	  }
	  return matchFound;

	}// End of TrackIdToMCTruth

	bool TrackMatchIdToMCTruth( Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth ){
	  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	  bool matchFound = false;
	  try {
	    art::Ptr<simb::MCTruth> temp_mctruth;
	    temp_mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
	    if (temp_mctruth == mctruth) {
	      matchFound = true;
	    }
	  } catch(...) {
	    std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
	    std::cout<<"Exception inside TrackMatchIdToMCTruth " << std::endl;;
	    matchFound = false;
	  }
	  return matchFound;

	}// End of TrackMatchIdToMCTruth



private:


};

} // end namespace analysis

#endif
