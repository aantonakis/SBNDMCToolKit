
set(
        MODULE_LIBRARIES
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                larsim::Simulation
                lardataobj::Simulation
                larsim::MCCheater_BackTrackerService_service
                larsim::MCCheater_ParticleInventoryService_service
                lardata::Utilities
                larevt::Filters
                lardataobj::RawData
                lardataobj::RecoBase
                lardata::RecoObjects
                larpandora::LArPandoraInterface
                nusimdata::SimulationBase
                sbndcode_CRT
                sbndcode_CRTData
                sbnobj::Common_CRT
                sbndcode_CRTUtils
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art::Persistency_Common
                art::Persistency_Provenance
                art::Utilities
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                ROOT::Geom
                ROOT::XMLIO
                ROOT::Gdml
                ROOT::Core
                sbndcode_RecoUtils
        )
cet_build_plugin(EtaFilter art::module SOURCE EtaFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

