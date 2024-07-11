import ROOT
import numpy as np
import sys
import os

print("Let's view something ...")

name = sys.argv[1]

print("ls file contents:")
inFile = ROOT.TFile.Open(name, "READ")
inFile.ls()

inFile.cd("ana")

inFile.ls()

tree = inFile.Get("ana/mc_nu_ttree")

tree.Print()


part_tree = inFile.Get("ana/mc_nu_particle_ttree")

part_tree.Print()


count = 0
for event in part_tree:
	EventID = getattr(event, "EventID")
	#nu_id = getattr(event, "NuID")
	pdg_vec = getattr(event, "genie_primaries_pdg")
	
	if count > 2:
		break

	count += 1

	print("EventID", EventID)
	#print("nu_id", nu_id)
	print("pdg_vec", pdg_vec)
	for num in pdg_vec:
		print("PDG ", num)


