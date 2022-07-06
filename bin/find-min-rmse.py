####!/usr/bin/env python3
#Author: Matthias Vogler

import numpy as np
import sys
import os

def usage():
	print("Usage: python3 find-min-rmse.py <BINDIR> <FRAGDIR> <WORKDIR> <NFRAGS> <NFIT> <MINCHGS> <MAXCHGS> <MINFRAGCHGS> <MAXFRAGCHGS>")


if len(sys.argv) < 10:
	usage()
	exit()

BINDIR = str(sys.argv[1])
FRAGDIR = str(sys.argv[2])
WORKDIR = str(sys.argv[3])
nfrag = int(sys.argv[4])
nfit = int(sys.argv[5])
minchgs = int(sys.argv[6])
maxchgs = int(sys.argv[7]) + 1
minfragchgs = int(sys.argv[8])
maxfragchgs = int(sys.argv[9])



def read_rmse(nfrag, nfit, minfragchgs, maxfragchgs):
	rmse_data = np.zeros([nfrag, maxfragchgs - minfragchgs + 1, nfit])
	for ifrag in range(0, nfrag):
		for ifit in range(0, nfit):
			for ichgs in range(0, maxfragchgs+1 - minfragchgs):
				xyzfile = FRAGDIR + "/frag" + str(ifrag+1) + "/fit" + str(ifit+1) + "/" + str(ichgs+minfragchgs) + "charges.xyz"
				#print(xyzfile)
				try:
					with open(xyzfile, 'r') as file:
						data = file.read().replace('\n', ' ').partition("RMSE")
						data = data[2].partition("kcal/mol")
						data = data[0].strip()
					rmse_data[ifrag, ichgs, ifit] = float(data)
				except FileNotFoundError as e:
					print("ERROR:" + xyzfile + " does not exist! \n Skipping...")
					rmse_data[ifrag, ichgs, ifit] = 1000.0
	return rmse_data


def combine_frags(nfrag, nfit, fit_vec, perm_vec, fragdir, BINDIR):
	base = "python3 " + BINDIR + "/combine-frags.py "
	string = ""
	for ifrag in range(0, nfrag):
		string += fragdir + "/frag" + str(ifrag+1) + "/fit" + str(int(fit_vec[ifrag] + 1)) + "/" + str(int(perm_vec[ifrag])) + "charges.xyz "
	return base + string


rmse_data = read_rmse(nfrag, nfit, minfragchgs, maxfragchgs)

#print(combine_frags(4, 3, [0, 2, 1, 0], [5, 4, 4, 4], fragdir))


for i in range(minchgs, maxchgs):
	os.system("python3 " + BINDIR + "/charge-permutations.py " + str(nfrag) + " " + str(i) + " " + str(minfragchgs) + " " + str(maxfragchgs) + " > permutations-" + str(i) + ".dat")
	if os.stat("permutations-" + str(i) + ".dat").st_size != 0:
		perm = np.loadtxt("permutations-" + str(i) + ".dat")
		nperm = perm.shape[0]
		best_rmse = 99999.0
		rmse_perm = np.zeros(nperm)
		fit_index = np.zeros([nfrag, nperm])
		for iperm in range(0, nperm):
			rmse_avg = 0.0
			for ifrag in range(0, nfrag):
				rmse = 99999.0
				for ifit in range(0, nfit):
					try:
						rmse_curr = rmse_data[ifrag, int(perm[iperm, ifrag]-minfragchgs), ifit]
					except IndexError as e:
						rmse_curr = rmse_data[ifrag, int(perm[iperm]-minfragchgs), ifit]
					if rmse_curr < rmse:
						rmse = rmse_curr
						fit_index[ifrag, iperm] = ifit
				rmse_avg += rmse
			rmse_perm[iperm] = rmse_avg/nfrag
		try:
			print("Found minimal solution with RMSE: "+str(np.min(rmse_perm)))
			os.system(combine_frags(nfrag, nfit, fit_index[:, np.argmin(rmse_perm)], perm[np.argmin(rmse_perm), :], FRAGDIR, BINDIR))
		except IndexError as e:
			print("Found minimal solution with RMSE: "+str(np.min(rmse_perm)))
			os.system(combine_frags(nfrag, nfit, fit_index[:, np.argmin(rmse_perm)], perm, FRAGDIR, BINDIR))		
		os.system("mv combined.xyz "+ WORKDIR + "/" + str(int(i)) + "-combined.xyz")
		os.system("rm permutations-" + str(i) + ".dat")
	else:
		print("WARNING: NO SOLUTION FOUND FOR " + str(i) + " charges!")
		os.system("rm permutations-" + str(i) + ".dat")




	

	
