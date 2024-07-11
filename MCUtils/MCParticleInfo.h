#ifndef MC_PARTICLE_INFO_H
#define MC_PARTICLE_INFO_H

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


namespace analysis {


class MCParticleInfo {

public:
        MCParticleInfo() {}
        ~MCParticleInfo() {}


	/*
 *	Currently, the purpose of this class is to grab both the Genie and Geant
 *	Particle information associated with a given MC neutrino and store that information 
 *	neutrino by neutrino. May want to rename this class to make that more apparent. We
 *	may want a similar class for cosmics for example and will want to easily distinguish 
 *	between the different trees
 * 	*/


	art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
	
	// Genie
	std::vector<int>   genie_nu_ids;
	std::vector<Int_t> genie_primaries_pdg;
	std::vector<float> genie_Eng;
	std::vector<float> genie_Px;
	std::vector<float> genie_Py;
	std::vector<float> genie_Pz;
	std::vector<float> genie_P;
	std::vector<Int_t> genie_status_code;
	std::vector<float> genie_mass;
	std::vector<Int_t> genie_trackID;
	std::vector<Int_t> genie_num_daughters;
	std::vector<Int_t> genie_mother;

        unsigned int fEventID;	
  	unsigned int fRun;
	unsigned int fSubRun;
	//unsigned int fNuID; // we will save the info neutrino by neutrino

	// Geant
	std::vector<int>  geant_nu_ids;
	std::vector<int>  geant_process_primary;
        std::vector<std::string> geant_processname;
        std::vector<std::string> geant_Endprocessname;
        std::vector<int> geant_Mother;
        std::vector<int> geant_TrackId;
        std::vector<int> geant_pdg;
        std::vector<int> geant_status;
        std::vector<float> geant_Eng;
        std::vector<float> geant_EndE;
        std::vector<float> geant_Mass;
        std::vector<float> geant_Px;
        std::vector<float> geant_Py;
        std::vector<float> geant_Pz;
        std::vector<float> geant_P;
        std::vector<float> geant_StartPointx;
        std::vector<float> geant_StartPointy;
        std::vector<float> geant_StartPointz;
        std::vector<float> geant_StartT;
        std::vector<float> geant_EndPointx;
        std::vector<float> geant_EndPointy;
        std::vector<float> geant_EndPointz;
        std::vector<float> geant_EndT;
        std::vector<float> geant_theta;
        std::vector<float> geant_phi;
        std::vector<float> geant_theta_xz;
        std::vector<float> geant_theta_yz;
        //geant_pathlen.push_back(plen);
        std::vector<int> geant_NumberDaughters;
        //geant_inTPCActive.push_back(int(isActive));
        //std::vector<float> geant_DepEnergy;
        //std::vector<float> geant_NumElectrons;

	// call this at the beginning of each neutrino interaction
	void clearMCParticleData() {

	  // Genie
	  genie_nu_ids.clear();
 	  genie_primaries_pdg.clear();
	  genie_Eng.clear();
	  genie_Px.clear();
	  genie_Py.clear();
	  genie_Pz.clear();
	  genie_P.clear();
	  genie_status_code.clear();
	  genie_mass.clear();
	  genie_trackID.clear();
	  genie_num_daughters.clear();
	  genie_mother.clear(); 

	  // Geant
	  geant_nu_ids.clear();
	  geant_process_primary.clear();
          geant_processname.clear();
	  geant_Endprocessname.clear();
          geant_Mother.clear();
          geant_TrackId.clear();
          geant_pdg.clear();
          geant_status.clear();
          geant_Eng.clear();
          geant_EndE.clear();
          geant_Mass.clear();
          geant_Px.clear();
          geant_Py.clear();
          geant_Pz.clear();
          geant_P.clear();
          geant_StartPointx.clear();
          geant_StartPointy.clear();
          geant_StartPointz.clear();
          geant_StartT.clear();
          geant_EndPointx.clear();
          geant_EndPointy.clear();
          geant_EndPointz.clear();
          geant_EndT.clear();
          geant_theta.clear();
          geant_phi.clear();
          geant_theta_xz.clear();
          geant_theta_yz.clear();
	  geant_NumberDaughters.clear();
	  //geant_DepEnergy.clear();
	  //geant_NumElectrons.clear();

	}

        void setMCParticleAddresses(TTree *tree) {

	  tree->Branch("EventID", &fEventID);
	  //tree->Branch("NuID", &fNuID);
	  tree->Branch("Run", &fRun);	  
	  tree->Branch("SubRun", &fSubRun);

	  // Genie
	  tree->Branch("genie_nu_ids", &genie_nu_ids);
	  tree->Branch("genie_primaries_pdg", &genie_primaries_pdg);
	  tree->Branch("genie_Eng", &genie_Eng);
	  tree->Branch("genie_Px", &genie_Px);
	  tree->Branch("genie_Py", &genie_Py);
	  tree->Branch("genie_Pz", &genie_Pz);
	  tree->Branch("genie_P", &genie_P);
	  tree->Branch("genie_status_code", &genie_status_code);
	  tree->Branch("genie_mass", &genie_mass);
	  tree->Branch("genie_trackID", &genie_trackID);
	  tree->Branch("genie_num_daughters", &genie_num_daughters);
	  tree->Branch("genie_mother", &genie_mother);
	
	  // Geant
	  tree->Branch("geant_nu_ids", &geant_nu_ids);
	  tree->Branch("geant_process_primary", &geant_process_primary);
	  tree->Branch("geant_pdg", &geant_pdg);
	  tree->Branch("geant_processname", &geant_processname);
	  tree->Branch("geant_Endprocessname", &geant_Endprocessname);
	  tree->Branch("geant_status", &geant_status);
	  tree->Branch("geant_Eng", &geant_Eng);
	  tree->Branch("geant_StartPointx", &geant_StartPointx);
	  tree->Branch("geant_StartPointy", &geant_StartPointy);
	  tree->Branch("geant_StartPointz", &geant_StartPointz);
	  tree->Branch("geant_EndPointx", &geant_EndPointx);
	  tree->Branch("geant_EndPointy", &geant_EndPointy);
	  tree->Branch("geant_EndPointz", &geant_EndPointz);
	  tree->Branch("geant_P", &geant_P);
	  tree->Branch("geant_Px", &geant_Px);
	  tree->Branch("geant_Py", &geant_Py);
	  tree->Branch("geant_Pz", &geant_Pz);
	  tree->Branch("geant_EndE", &geant_EndE);
	  tree->Branch("geant_EndT", &geant_EndT);
	  tree->Branch("geant_NumberDaughters", &geant_NumberDaughters);	
	  //tree->Branch("geant_DepEnergy", &geant_DepEnergy);
	  tree->Branch("geant_TrackId", &geant_TrackId);
	  tree->Branch("geant_Mass", &geant_Mass);
	  tree->Branch("geant_Mother", &geant_Mother);
	  tree->Branch("geant_StartT", &geant_StartT);
	  tree->Branch("geant_theta", &geant_theta);
	  tree->Branch("geant_theta_xz", &geant_theta_xz);
	  tree->Branch("geant_theta_yz", &geant_theta_yz);
	  tree->Branch("geant_phi", &geant_phi);
	  //tree->Branch("geant_NumElectrons", &geant_NumElectrons);


	}



        //void fillMCParticleTree(TTree *tree, art::Event const& e, art::Ptr<simb::MCTruth>& mctruth) {
        void fillMCParticleTree(TTree *tree, art::Event const& e, std::string nugenlabel) {
	  	 
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();
	  //fNuID = nu_id;
	  
	  // * MC truth information --> How we will get Neutrino Info
  	  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  	  std::vector<art::Ptr<simb::MCTruth> > mclist;
  	  if (e.getByLabel(nugenlabel, mctruthListHandle))
  	    art::fill_ptr_vector(mclist, mctruthListHandle);
	
	  unsigned int nu_id = 0;

	  // Loop over the MC truths and find the neutrinos
	  for (unsigned int i_mctruth = 0; i_mctruth < mclist.size(); i_mctruth++){
            // fetch a mctruth
            art::Ptr<simb::MCTruth> mctruth = mclist[i_mctruth];
            //Check if it's a neutrino
            if (!mctruth->NeutrinoSet()) continue;


	    // Loop over the Genie particles for a given neutrino
	    size_t StoreParticles = mctruth->NParticles();
	    //std::cout << "Number of Genie Particles " << StoreParticles << std::endl;
	    for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
              const simb::MCParticle& part(mctruth->GetParticle(iPart));
	      //std::cout << "primary pdg " << part.PdgCode() << std::endl;
              genie_nu_ids.push_back(nu_id);
	      genie_primaries_pdg.push_back(part.PdgCode());
              genie_Eng.push_back(part.E());
              genie_Px.push_back(part.Px());
              genie_Py.push_back(part.Py());
              genie_Pz.push_back(part.Pz());
              genie_P.push_back(part.P());
              genie_status_code.push_back(part.StatusCode());
              genie_mass.push_back(part.Mass());
              genie_trackID.push_back(part.TrackId());
              genie_num_daughters.push_back(part.NumberDaughters());
              genie_mother.push_back(part.Mother());

            } // end Genie particle loop


	    const sim::ParticleList& plist = pi_serv->ParticleList();

	    std::string pri("primary");
      	    int primary=0;
      	    //int active = 0;
      	    int geant_particle=0;
      	    sim::ParticleList::const_iterator itPart = plist.begin(),
            pend = plist.end(); // iterator to pairs (track id, particle)

	    // Loop over geant particles and find if they belong to the nu interaction

      	    for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
              const simb::MCParticle* pPart = (itPart++)->second;
              if (!pPart) {
                throw art::Exception(art::errors::LogicError)
                << "GEANT particle #" << iPart << " returned a null pointer";
              }
	      // Grab the mctruth associated with this particle
	      const art::Ptr<simb::MCTruth> matched_mctruth = pi_serv->ParticleToMCTruth_P(pPart);
	      // check if it's from our neutrino
	      if (!(matched_mctruth == mctruth))
	        continue;
 
	      ++geant_particle;

	      bool isPrimary = pPart->Process() == pri;
              if (isPrimary) ++primary;
              int TrackID = pPart->TrackId();
	      geant_nu_ids.push_back(nu_id);
	      geant_process_primary.push_back(int(isPrimary));
              geant_processname.push_back(pPart->Process());
              geant_Endprocessname.push_back(pPart->EndProcess());
              geant_Mother.push_back(pPart->Mother());
              geant_TrackId.push_back(TrackID);
              geant_pdg.push_back(pPart->PdgCode());
              geant_status.push_back(pPart->StatusCode());
              geant_Eng.push_back(pPart->E());
              geant_EndE.push_back(pPart->EndE());
              geant_Mass.push_back(pPart->Mass());
              geant_Px.push_back(pPart->Px());
              geant_Py.push_back(pPart->Py());
              geant_Pz.push_back(pPart->Pz());
              geant_P.push_back(pPart->Momentum().Vect().Mag());
              geant_StartPointx.push_back(pPart->Vx());
              geant_StartPointy.push_back(pPart->Vy());
              geant_StartPointz.push_back(pPart->Vz());
              geant_StartT.push_back(pPart->T());
              geant_EndPointx.push_back(pPart->EndPosition()[0]);
              geant_EndPointy.push_back(pPart->EndPosition()[1]);
              geant_EndPointz.push_back(pPart->EndPosition()[2]);
              geant_EndT.push_back(pPart->EndT());
              geant_theta.push_back(pPart->Momentum().Theta());
              geant_phi.push_back(pPart->Momentum().Phi());
              geant_theta_xz.push_back(std::atan2(pPart->Px(), pPart->Pz()));
              geant_theta_yz.push_back(std::atan2(pPart->Py(), pPart->Pz()));
              //geant_pathlen.push_back(plen);
              geant_NumberDaughters.push_back(pPart->NumberDaughters());
              //geant_inTPCActive.push_back(int(isActive));
              //geant_DepEnergy.push_back(totalE_particle);
              //geant_NumElectrons.push_back(totalne_particle);

	    } // end geant particle loop
	    nu_id += 1;
	  } // end of loop over mctruths  
          tree->Fill();

	} // end of fill function


private:


};

} // end namespace analysis

#endif


