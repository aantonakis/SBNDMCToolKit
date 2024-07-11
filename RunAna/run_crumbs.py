import os
import sys

file_list = sys.argv[1]
out_name = sys.argv[2]

os.system('lar -c crumbs_sbnd.fcl -S '+file_list+" -o /pnfs/sbnd/scratch/users/aantonak/EtaSim/"+out_name)


