import os

datasetLst = [line.rstrip('\n') for line in open('NCTC_ds.txt')]
for i in range(len(datasetLst)):
	if len(datasetLst[i])==0:
		break
	ds = datasetLst[i]
	fldr = ds + "_filtered"

	os.system("mkdir -p {}".format(fldr))
	os.system("cp filtered_fasta/{}_reads.fasta {}/reads.fasta".format(ds,fldr))
	os.system("cp groundTruths/{}_daligner_ground_truth.txt {}/daligner_ground_truth.txt".format(ds,fldr))
