#ifndef RECO_SHOWER_INFO_H
#define RECO_SHOWER_INFO_H

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


class RecoShowerInfo {

public:
        RecoShowerInfo() {}
        ~RecoShowerInfo() {}

        unsigned int fEventID;
  	unsigned int fRun;
	unsigned int fSubRun;
        // Attributes of Reconstructed Showers
	Short_t  nshowers;                  // Number of reconstructed showers
        std::vector<Int_t>  fshwID;        // shower ID
	std::vector<Int_t>  fshwBestPlane;
        //std::vector<std::vector<double>> fshwdEdx;        
        //std::vector<std::vector<double>> fshwdEdxErr; 
	std::vector<double> fshwDirX;       
	std::vector<double> fshwDirY;       
	std::vector<double> fshwDirZ;       
	std::vector<double> fshwDirXErr;       
	std::vector<double> fshwDirYErr;       
	std::vector<double> fshwDirZErr;       
	std::vector<std::vector<double>> fshwEnergy;
	std::vector<std::vector<double>> fshwEnergyErr;
	std::vector<double> fshwLength;
	std::vector<double> fshwOpenAngle;
	std::vector<double> fshwStartX;       
	std::vector<double> fshwStartY;       
	std::vector<double> fshwStartZ;       
	std::vector<double> fshwStartXErr;       
	std::vector<double> fshwStartYErr;       
	std::vector<double> fshwStartZErr;       
	std::vector<Int_t>  fshwPFPID; // ID of the corresponding PFP
	std::vector<Int_t>  fshwIsPrimary; // Is it Primary according to Reco
	std::vector<Int_t>  fshwNDaughters; // How many daughter PFPs
	std::vector<Int_t>  fshwParentPFPID; // The parent PFP ID

        typedef std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle> > shwPfpMap;
	typedef std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>>::const_iterator shwPfpMapIt;
	shwPfpMap showerPFParticleMap;

        // call this at the beginning of each LArSoft Event
	void clearRecoShowerData() {
          fshwID.clear();
	  fshwBestPlane.clear();
	  fshwEnergy.clear();
	  fshwEnergyErr.clear();
	  fshwLength.clear();
	  fshwOpenAngle.clear();
	  fshwDirX.clear(); fshwDirY.clear(); fshwDirZ.clear();
	  fshwDirXErr.clear(); fshwDirYErr.clear(); fshwDirZErr.clear();
	  fshwStartX.clear(); fshwStartY.clear(); fshwStartZ.clear();
	  fshwStartXErr.clear(); fshwStartYErr.clear(); fshwStartZErr.clear();
	  fshwPFPID.clear();
	  fshwIsPrimary.clear();
	  fshwNDaughters.clear();
	  fshwParentPFPID.clear();
	  showerPFParticleMap.clear();

        } // end of clearRecoShowerData

        void setRecoShowerAddresses(TTree *tree,  std::string showerlabel) {
          std::string ShowerLabel = showerlabel;
 	  std::string BranchName = "shwID_"+ShowerLabel; 
  	  std::cout << "BranchName " << BranchName << std::endl;
	  tree->Branch("EventID", &fEventID);
	  tree->Branch("Run", &fRun);	  
	  tree->Branch("SubRun", &fSubRun);
	  tree->Branch("shwID", &fshwID);
	  tree->Branch("shwBestPlane", &fshwBestPlane);
	  tree->Branch("shwEnergy", &fshwEnergy);
	  tree->Branch("shwEnergyErr", &fshwEnergyErr);
	  tree->Branch("shwLength", &fshwLength);
	  tree->Branch("shwOpenAngle", &fshwOpenAngle);
	  tree->Branch("shwDirX", &fshwDirX);
	  tree->Branch("shwDirY", &fshwDirY);
	  tree->Branch("shwDirZ", &fshwDirZ);
	  tree->Branch("shwDirXErr", &fshwDirXErr);
	  tree->Branch("shwDirYErr", &fshwDirYErr);
	  tree->Branch("shwDirZErr", &fshwDirZErr);
	  tree->Branch("shwStartX", &fshwStartX);
	  tree->Branch("shwStartY", &fshwStartY);
	  tree->Branch("shwStartZ", &fshwStartZ);
	  tree->Branch("shwStartXErr", &fshwStartXErr);
	  tree->Branch("shwStartYErr", &fshwStartYErr);
	  tree->Branch("shwStartZErr", &fshwStartZErr);
	  tree->Branch("shwPFPID", &fshwPFPID);
	  tree->Branch("shwIsPrimary", &fshwIsPrimary);
	  tree->Branch("shwNDaughters", &fshwNDaughters);
	  tree->Branch("shwParentPFPID", &fshwParentPFPID);
        
        } // end of setRecoShowerAddresses


        void makeShwPFPMap(art::Event const& e, std::string showerlabel, art::Handle<std::vector<recob::PFParticle>> pfpHandle) {
      
	  art::FindManyP<recob::Shower> fshw(pfpHandle, e, showerlabel);
          
	  for(unsigned int i = 0; i < pfpHandle->size(); ++i) {
            art::Ptr< recob::PFParticle > pfp( pfpHandle, i );

            // Get shower association
            std::vector< art::Ptr<recob::Shower> > shwAssn = fshw.at(pfp->Self());
                 
            if(shwAssn.size()  > 1){
              std::cout << "PFParticle has " << shwAssn.size() << " associated showers, should only have 1 or 0 ";
              continue;
            }
            if(shwAssn.size() == 0) continue;
            showerPFParticleMap.emplace(shwAssn[0],pfp);

          }
         
        } // end of makeShwPFPMap

        void fillRecoShowerTree(TTree *tree, art::Event const& e, std::string showerlabel, std::string pfplabel) {
        
	  // Want to record what event we are in
	  fEventID = e.id().event();
  	  fRun = e.run();
  	  fSubRun = e.subRun();

	  // Load showers from Pandora or a different producer
  	  art::Handle<std::vector<recob::Shower>> showerHandle;
          std::vector<art::Ptr<recob::Shower>> showerVec;
  	  if (e.getByLabel(showerlabel, showerHandle))
            art::fill_ptr_vector(showerVec, showerHandle);
          
	  if (showerVec.empty())
            return;

	  std::cout << "The PFP label being matched to showers " << pfplabel << std::endl;
	  //art::FindManyP<recob::PFParticle> fm_shwPFP(showerVec, e, pfplabel);

          // Loop over all of the reco showers in this event
    	  for (const art::Ptr<recob::Shower> &shw : showerVec) {
  	    shwPfpMapIt it;
            it = showerPFParticleMap.find(shw);
            if(it == showerPFParticleMap.end()) continue; // check if it's in the showerPFP map
            art::Ptr<recob::PFParticle> pfpAssn = it->second;

	    // First, let's match the pfp to this shower
	    //std::vector<art::Ptr<recob::PFParticle>> pfpAssn = fm_shwPFP.at(shw.key());
            //if(pfpAssn.size() == 1) {
	    fshwPFPID.push_back(pfpAssn->Self());
	    fshwIsPrimary.push_back(pfpAssn->IsPrimary());
	    fshwNDaughters.push_back(pfpAssn->NumDaughters());
	    fshwParentPFPID.push_back(pfpAssn->Parent());
	    //}
            fshwID.push_back(shw->ID()); 
	    fshwBestPlane.push_back(shw->best_plane()); 
	    fshwEnergy.push_back(shw->Energy()); 
	    fshwEnergyErr.push_back(shw->EnergyErr());
	    
	    if (shw->has_open_angle()) 
	      fshwOpenAngle.push_back(shw->OpenAngle()); 
	    else
	      fshwOpenAngle.push_back(-1.);
 
	    if (shw->has_length()) 
	      fshwLength.push_back(shw->Length()); 
	    else
	      fshwLength.push_back(-1.); 
	 
	    fshwDirX.push_back(shw->Direction().X());	
	    fshwDirY.push_back(shw->Direction().Y());	
	    fshwDirZ.push_back(shw->Direction().Z());	
	    fshwDirXErr.push_back(shw->DirectionErr().X());	
	    fshwDirYErr.push_back(shw->DirectionErr().Y());	
	    fshwDirZErr.push_back(shw->DirectionErr().Z());	
	    fshwStartX.push_back(shw->ShowerStart().X());
	    fshwStartY.push_back(shw->ShowerStart().Y());
	    fshwStartZ.push_back(shw->ShowerStart().Z());
	    fshwStartXErr.push_back(shw->ShowerStartErr().X());
	    fshwStartYErr.push_back(shw->ShowerStartErr().Y());
	    fshwStartZErr.push_back(shw->ShowerStartErr().Z());

          } // end of loop over showers

	  // TODO Note: we could put the Fill statement in the analyzer module if we decide to
	  // incorporate more info outside of this simple loop

          tree->Fill();
	 
        } // end of fillRecoShowerTree


private:
  // might not need this


};

} // end namespace analysis

#endif






