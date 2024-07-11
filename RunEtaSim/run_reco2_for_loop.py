import sys
import os

file_list = sys.argv[1]

out_dir = "/pnfs/sbnd/scratch/users/aantonak/EtaSim/Reco2_Batch2/"


count = 0
with open(file_list, 'r') as file:
	for line in file:
	# Process each line
		reco1 = line.strip() 

		os.system("lar -c reco2_sce.fcl -s " + reco1 + " -o " + out_dir+"out_reco2_etasim_tpc_"+str(count)+".root")


		


