import os

os.system("lar -c prodgenie_nu_spill_tpc_sbnd_eta.fcl -n 6000 -o out_gen_test_6000.root")
os.system("lar -c g4_sce_lite.fcl -s out_gen_test_6000.root -o out_g4_test_6000.root")
os.system("lar -c detsim_sce_lite.fcl -s out_g4_test_6000.root -o out_detsim_test_6000.root")
os.system("lar -c reco1_sce.fcl -s out_detsim_test_6000.root -o out_reco1_test_6000.root")
os.system("lar -c reco2_sce.fcl -s out_reco1_test_6000.root -o out_reco2_test_6000.root")

