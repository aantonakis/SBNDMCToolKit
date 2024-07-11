import numpy as np
import os
from glob import iglob
import sys

batch = input("Give Batch:")



data_dir = "/pnfs/sbnd/scratch/users/aantonak/EtaSim/"+batch+"/**/*"

file_list = [f for f in iglob(data_dir, recursive=True) if os.path.isfile(f)]	


print("Extracting All Files ...", "\n")
for f in file_list:
	print(f)

print("\n", "Extracting Simulation ROOT files")

sim_files = []
for f in file_list:
	if "prodgenie" in f and "hists" not in f and "json" not in f:
		sim_files.append(f)
		print(f)

print('\n')
print("Number of Files for next Step:", len(sim_files))
print('\n')


print("Writing to a textfile")
with open('reco1_'+batch+'_'+'nf_'+str(len(sim_files))+'_list.txt', 'w') as F:
	for f in sim_files:
		F.write(f)
		F.write('\n')






