# Generates minHashes for a dataset and stores the corresponding minHash Matrix in dataset/minHashes/minHashArr.txt

from Bio.SeqIO.FastaIO import SimpleFastaParser
import mmh3
import itertools
import multiprocessing as mp
import numpy as np
import os, pickle
import argparse
from tqdm import tqdm


def runSim(args):
	idx,read = args[0]
	dataset, numHashes = args[1:]

	symLength = 7  # k in k-mer

	## create list of k-mers
	iterLst = []
	for j in range(len(read) - symLength):
		iterLst.append(read[j:j + symLength])

	hashOut = []
	hashOutArgMin = []
	for j in range(numHashes):
		hashLst = [mmh3.hash(sym, j, signed=False) for sym in iterLst]
		minVal = min(hashLst)
		hashOut.append(minVal)
		aMin = np.argmin(hashLst)
		hashOutArgMin += [(aMin,iterLst[aMin])]

	minHashes = np.array(hashOut)
	minHashes.dump(dataset+"minHashes/minHashes_{}.txt".format(idx))
	pickle.dump( hashOutArgMin, open( dataset+"minHashes/argminHashes_{}.txt".format(idx), "wb" ) )


ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. NCTC5047_filtered")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--num_hashes",  help="Number of hashes to compute", type=int,default=1000)

args = ap.parse_args()
num_jobs   = args.num_jobs
dataset    = args.dataset        
numHashes  = args.num_hashes

if dataset[-1] != '/':
	dataset += '/'


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])
n= len(reads_lst)

if not os.path.isdir(dataset+"/minHashes"):
	print("Creating minHashes folder")
	os.mkdir(dataset+"/minHashes")

print("finished loading reads \n parallelizing minHash computations")

pool      = mp.Pool(processes=num_jobs)
arg_tuple =  itertools.product(enumerate(reads_lst), [dataset], [numHashes])
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=n):
	pass

## save all minHashes in one array
minHashArr = np.zeros((n,numHashes))
for i in range(n):
    minHashArr[i] = np.load(dataset + "minHashes/minHashes_{}.txt".format(i))
print("writing minHashArr")
np.savetxt(dataset + "minHashes/minHashArr.txt",minHashArr)

## generate dataset specific random calibration reads
## generate 1 per thread, can increase this to have more random reads
cmd = "python generateRandReads.py --dataset {} --num_jobs {} --num_hashes {} --num_reads {}".format(
	dataset,num_jobs,numHashes,num_jobs)
print("running: ",cmd)
os.system(cmd)

## precompute all jaccard similarities for comparison
cmd = "python computeJSim.py --dataset {} --num_jobs {}".format(
	dataset,num_jobs)
print("running: ",cmd)
os.system(cmd)
