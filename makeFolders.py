import os
import argparse

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--datasets",  help="Text file with folder to dataset on each line", type = str, default = "NCTC_ds.txt")
args = ap.parse_args()
datasets = args.datasets


datasetLst = [line.rstrip('\n') for line in open(datasets)]
for i in range(len(datasetLst)):
	if len(datasetLst[i])==0:
		break
	ds = datasetLst[i]
	fldr = ds + "_filtered"

	os.system("mkdir -p {}".format(fldr))
	os.system("cp filtered_fasta/{}_reads.fasta {}/reads.fasta".format(ds,fldr))
	os.system("cp groundTruths/{}_daligner_ground_truth.txt {}/daligner_ground_truth.txt".format(ds,fldr))
