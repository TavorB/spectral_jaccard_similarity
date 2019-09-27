from Bio.SeqIO.FastaIO import SimpleFastaParser
import itertools
import multiprocessing as mp
import numpy as np
import os
import argparse
import shutil
from tqdm import tqdm

def runSim(args):

	## avoid one processes starting multiple threads
	os.environ["MKL_NUM_THREADS"] = "1" 
	os.environ["NUMEXPR_NUM_THREADS"] = "1" 
	os.environ["OMP_NUM_THREADS"] = "1" 

	ref_read = args[0]
	dataset,numRandReads = args[1]


	## load precomputed minHashes
	minHashArr = np.loadtxt(dataset+"minHashes/minHashArr.txt")
	randMinHashArr = np.zeros((numRandReads,1000))
	for i in range(numRandReads):
		randMinHashArr[i] = np.load(dataset+"randReads/randMinHashes_{}.txt".format(i))
	minHashArrExtended = np.vstack((minHashArr[:,:1000],randMinHashArr))

	## minHash collision matrix
	empiricalMatrix = (minHashArrExtended == minHashArrExtended[ref_read])
	empiricalMatrix = np.delete(empiricalMatrix,ref_read,axis=0)


	## compute SVd
	U,s,VT = np.linalg.svd(empiricalMatrix - np.ones(empiricalMatrix.shape))

	u = U[:,0]
	pHatSVD = 1-np.abs(u[:n-1])/np.abs(np.median(u[n-1:]))
	
	qTemp = VT[0,:]
	qs = 1-np.abs(qTemp)/np.max(np.abs(qTemp))


	## save results
	np.savetxt(dataset+"/SVD/pi_refread_{}.txt".format(ref_read),pHatSVD)
	np.savetxt(dataset+"/SVD/qj_refread_{}.txt".format(ref_read),qs)

	np.savetxt(dataset+"/SVD/raw_pi_refread_{}.txt".format(ref_read),u)
	np.savetxt(dataset+"/SVD/raw_qj_refread_{}.txt".format(ref_read),qTemp)


ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--num_rand", help="Number of random calibration reads", type = int, default=5)

args = ap.parse_args()

num_jobs   = args.num_jobs
dataset    = args.dataset        
numRandReads = args.num_rand

if dataset[-1] != '/':
	dataset += '/'

## delete all old SVDs
if os.path.isdir(dataset+"/SVD"):
	print("Deleting all old SVDs")
	shutil.rmtree(dataset+"/SVD")
print("Creating SVD folder")
os.mkdir(dataset+"/SVD")


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])
n = len(reads_lst)

print("Parallelizing SVDs")
pool      = mp.Pool(processes=num_jobs)
arg_tuple =  itertools.product(range(n), [[dataset,numRandReads]])
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=n):
	pass