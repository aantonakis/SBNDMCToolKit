////////////////////////////////////////////////////////////////////////
// Class:       MCScraper
// Plugin Type: analyzer (Unknown Unknown)
// File:        MCScraper_module.cc
//
// Generated at Tue Jan 23 01:17:55 2024 by Alexander Antonakis using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "canvas/Persistency/Common/FindManyP.h"

// Additional LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
//#include "sbnobj/Common/Reco/CNNScore.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <utility>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> // std::mem_fn()
#include <typeinfo>
#include <cmath>

#include "TTimeStamp.h"

#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <bitset>

#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"


// includes inspired from uboone code

#include "TDatabasePDG.h"
#include "TParticlePDG.h"


// Custom Classes
#include "MCUtils/TruthMatcher.h"
#include "MCUtils/GeoHelper.h"
#include "MCUtils/MCNeutrinoInfo.h"
#include "MCUtils/MCParticleInfo.h"
#include "MCUtils/RecoShowerInfo.h"
#include "MCUtils/RecoTrackInfo.h"
#include "MCUtils/PFParticleInfo.h"
#include "MCUtils/SliceInfo.h"
#include "MCUtils/CrumbsInfo.h"

namespace analysis {
  class MCScraper;

  /// information from the subrun
  struct SubRunData_t {
    SubRunData_t() { Clear(); }
    void Clear() { pot = -99999.; }
    Double_t pot; //protons on target
    Double_t goodpot; // good pots
  }; // struct SubRunData_t

}


class analysis::MCScraper : public art::EDAnalyzer {
public:
  explicit MCScraper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCScraper(MCScraper const&) = delete;
  MCScraper(MCScraper&&) = delete;
  MCScraper& operator=(MCScraper const&) = delete;
  MCScraper& operator=(MCScraper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(const art::SubRun& sr);
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  //art::InputTag potSummaryTag_;
  // Declare member data here.

  // Create a SubRun TTree
  TTree *fSubRunTree;  

  // Create an output TTree
  TTree *fTree;

  SubRunData_t SubRunData;

  Double_t totPOT = 0.;
  Double_t totGoodPOT = 0.;
 
  TTree *fpfpTree; 
  TTree *fSliceTree;
  TTree *fShwTree;
  TTree *fTrkTree;
  TTree *fnuTree;
  TTree *fnuPartTree;
  TTree *fCrumbsTree;


  PFParticleInfo pfpInfo;
  SliceInfo sliceInfo;
  RecoShowerInfo recoShwInfo; // custom object to handle reco showers
  RecoTrackInfo recoTrkInfo;
  MCNeutrinoInfo mcNuInfo; // custom object to handle our favorite particle -> nu
  MCParticleInfo mcNuPartInfo;
  CrumbsInfo crumbsInfo;

  // Tree variables
  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  double fEventTime;


  // Define input labels --> producers
  const std::string fPFParticleLabel; // Pandora for instance
  const std::string fTrackLabel;
  const std::string fShowerLabel; 
  const std::string fCalorimetryLabel; 
  const std::string fPOTModuleLabel;
  const std::string fVertexModuleLabel;
  const std::string fGenieGenLabel;
  const std::string fHitModuleLabel;
  const std::string fCosmicTaggerLabel;
  const std::string fOpT0FinderLabel;
  const std::string fCrumbsResultLabel;

};


analysis::MCScraper::MCScraper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")), // Put all of the producer labels for DataProducts here
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fPOTModuleLabel(p.get<std::string>("POTModuleLabel")),
  fVertexModuleLabel(p.get<std::string>("VertexModuleLabel")),
  fGenieGenLabel(p.get<std::string>("GenieGenLabel")),
  fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
  fCosmicTaggerLabel(p.get<std::string>("CosmicTaggerLabel")),
  fOpT0FinderLabel(p.get<std::string>("OpT0FinderLabel")),
  fCrumbsResultLabel(p.get<std::string>("CrumbsResultLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


void analysis::MCScraper::beginSubRun(const art::SubRun& sr)
{
  std::cout << "In the beginSubRun Function" << std::endl;
  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);

  if(sr.getByLabel(fPOTModuleLabel,potListHandle)) {
    SubRunData.pot=potListHandle->totpot;
    SubRunData.goodpot=potListHandle->totgoodpot;
    totPOT += potListHandle->totpot;
    totGoodPOT += potListHandle->totgoodpot;
    std::cout << "POT in SubRun " << potListHandle->totpot << std::endl;
    std::cout << "Good POT in SubRun " << potListHandle->totgoodpot << std::endl;
  }
  else {
    SubRunData.pot=0.;
  }
} // End of beginSubRun function



void analysis::MCScraper::analyze(art::Event const& e)
{
  // Grab the LArSoft Event ID for the current event
  fEventID = e.id().event();
  std::cout << "Processing Event: " << fEventID << std::endl;  

  fRun = e.run();
  fSubRun = e.subRun();

  std::cout << "Processing SubRun " << fSubRun << std::endl;

  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fEventTime = tts.AsDouble();


  //copied from MergeDataPaddles.cxx
  //  art::Handle< raw::BeamInfo > beam;
  //    if (evt.getByLabel("beam",beam)){
  //        fData->beamtime = (double)beam->get_t_ms();
  //            fData->beamtime/=1000.; //in second
  //              }
  //

 
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventoryService;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  art::Handle<std::vector<anab::T0>> T0Handle;
  std::vector<art::Ptr<anab::T0>> T0Vec;
  if (e.getByLabel(fPFParticleLabel, T0Handle))
    art::fill_ptr_vector(T0Vec, T0Handle);


  // Load the PFParticles from Pandora
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  // Check if there are any reconstructed particles in this event
  if (pfpVec.empty())
    return;

  // Custom class functions
  pfpInfo.clearPFParticleData();
  pfpInfo.fillPFParticleTree(fpfpTree, e, fPFParticleLabel, pfpHandle, fVertexModuleLabel, fOpT0FinderLabel);

  // Not Efficient, but convenient
  sliceInfo.clearSliceData();
  sliceInfo.fillSliceTree(fSliceTree, e, fPFParticleLabel, pfpHandle, fVertexModuleLabel, fOpT0FinderLabel);

  // custom class function
  recoShwInfo.clearRecoShowerData();
  // Make a map between showers and PFPs
  recoShwInfo.makeShwPFPMap(e, fShowerLabel, pfpHandle);
  // Fill the reco shower tree
  recoShwInfo.fillRecoShowerTree(fShwTree, e, fShowerLabel, fPFParticleLabel);

  recoTrkInfo.clearRecoTrackData();
  recoTrkInfo.makeTrkPFPMap(e, fTrackLabel, pfpHandle);
  recoTrkInfo.fillRecoTrackTree(fTrkTree, e, fTrackLabel, fPFParticleLabel, fCalorimetryLabel, fCosmicTaggerLabel);

  mcNuInfo.clearMCNeutrinoData(); 
  mcNuInfo.fillMCNeutrinoTree(fnuTree, e, fGenieGenLabel, fHitModuleLabel, fPFParticleLabel, clockData);

  mcNuPartInfo.clearMCParticleData();
  mcNuPartInfo.fillMCParticleTree(fnuPartTree, e, fGenieGenLabel);

  crumbsInfo.clearCrumbsData();
  crumbsInfo.fillCrumbsTree(fCrumbsTree, e, fPFParticleLabel, fCrumbsResultLabel);

  fTree->Fill();

} // end of analyze event

void analysis::MCScraper::beginJob()
{
  
  // Implementation of optional member function here.
  std::cout << "/-------- Begin Job Step -------- /" << std::endl; 
  
  std::cout << std::endl;
  std::cout << "Clearing the SubRun Info at the Begin Job Step ???" << std::endl;
  SubRunData.Clear();

  art::ServiceHandle<art::TFileService> tfs;
 

  fSubRunTree = tfs->make<TTree>("subrun_ttree", "subrun_ttree");
  fSubRunTree->Branch("totPOT", &totPOT);
  fSubRunTree->Branch("totGoodPOT", &totGoodPOT);
 
  fpfpTree = tfs->make<TTree>("pfp_ttree", "pfp_ttree");
  pfpInfo.setPFParticleAddresses(fpfpTree, fPFParticleLabel);
  
  fSliceTree = tfs->make<TTree>("slice_ttree", "slice_ttree");
  sliceInfo.setSliceAddresses(fSliceTree, fPFParticleLabel);

  fShwTree = tfs->make<TTree>("reco_shower_ttree", "reco_shower_ttree"); 
  recoShwInfo.setRecoShowerAddresses(fShwTree,  fShowerLabel);

  fTrkTree = tfs->make<TTree>("reco_track_ttree", "reco_track_ttree");
  recoTrkInfo.setRecoTrackAddresses(fTrkTree, fTrackLabel);  

  fnuTree = tfs->make<TTree>("mc_nu_ttree", "mc_nu_ttree");
  mcNuInfo.setMCNeutrinoAddresses(fnuTree);
 
  fnuPartTree = tfs->make<TTree>("mc_nu_particle_ttree", "mc_nu_particle_ttree");
  mcNuPartInfo.setMCParticleAddresses(fnuPartTree);

  fCrumbsTree = tfs->make<TTree>("crumbs_ttree", "crumbs_ttree");
  crumbsInfo.setCrumbsAddresses(fCrumbsTree);


  fTree = tfs->make<TTree>("tree", "output ttree");
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("Run", &fRun);
  fTree->Branch("SubRun", &fSubRun);
  fTree->Branch("EventTime", &fEventTime);
  
}

void analysis::MCScraper::endJob()
{
  // Implementation of optional member function here.

  // Fill the total POT Info at the end  
  fSubRunTree->Fill();
}

DEFINE_ART_MODULE(analysis::MCScraper)
