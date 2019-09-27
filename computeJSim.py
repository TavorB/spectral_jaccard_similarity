import argparse, os
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm,trange
import multiprocessing as mp
import shutil, itertools
from tqdm import tqdm


def generateSymSets(reads_lst,symLength):
	symSets = [0]*n
	for i,read in enumerate(reads_lst):
		lst = set()
		for j in range(len(read)-symLength):
			lst.add(read[j:j+symLength])
		symSets[i] = lst
	return symSets

def runSim(args):
	symLength=7 # k in k-mer
	ref_read,reads_lst = args

	symSets = generateSymSets(reads_lst,symLength)

	val = [len(symSets[ref_read].intersection(symSets[j]))/(len(symSets[ref_read].union(symSets[j]))) for j in range(ref_read+1,n)]
	np.savetxt(dataset+"/JSim_temp/JSim_{}.txt".format(ref_read),np.array(val))


ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )

args = ap.parse_args()
dataset    = args.dataset        
num_jobs   = args.num_jobs

if dataset[-1] != '/':
	dataset += '/'


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])
n = len(reads_lst)


if not os.path.isdir(dataset+"/JSim_temp"):
	os.mkdir(dataset+"/JSim_temp")

pool      = mp.Pool(processes=num_jobs)
arg_tuple =  itertools.product(range(n-1), [reads_lst])
print("Parallelizing JSim computation: ")
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=n-1):
	pass

## reconstruct JS array from individual files
JSims = np.eye(n)
print("Reconstructing JSim arr")
for ref_read1 in trange(n-1):
	val = np.loadtxt(dataset+"/JSim_temp/JSim_{}.txt".format(ref_read1))
	JSims[ref_read1,ref_read1+1:n] = val
	JSims[ref_read1+1:n,ref_read1] = val

if not os.path.isdir(dataset+"/minHashes"):
	print("Creating minHashes folder")
	os.mkdir(dataset+"/minHashes")
np.savetxt(dataset+"/minHashes/JSims.txt",JSims)

## delete all intermediate JS files
shutil.rmtree(dataset+"/JSim_temp")
print("deleting JSim_temp folder")
