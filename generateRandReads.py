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
	i = args[0]
	dataset, numHashes = args[1:]
	symLength = 7  # k in k-mer

	hashOut = []
	np.random.seed(i)
	rand_kmers = np.random.choice(list(kmer_dict.keys()),size = avgReadLen-symLength+1, replace=True, p=ps)
	
	for j in range(numHashes):
		hashLst = [mmh3.hash(sym, j, signed=False) for sym in rand_kmers]
		hashOut.append(min(hashLst))

	minHashes = np.array(hashOut)
	minHashes.dump(dataset+"randReads/randMinHashes_{}.txt".format(i))

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--num_reads",  help="Number of desired random calibration reads", type=int, default=5 )
ap.add_argument("--num_hashes",  help="Number of hashes to compute", type=int,default=1000)


args = ap.parse_args()

dataset    = args.dataset        
numReads   = args.num_reads
numHashes  = args.num_hashes
num_jobs   = min(args.num_jobs,numReads)

if dataset[-1] != '/':
	dataset += '/'


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])

n = len(reads_lst)
avgReadLen = int(np.mean([len(x) for x in reads_lst]))

if not os.path.isdir(dataset+"/randReads"):
	print("Creating randReads folder")
	os.mkdir(dataset+"/randReads")


## generate kmer frequencies for selecting kmers from dataset dist
print("Generating kmer frequencies")
kmer_dict = {}
k = 7
for ref_read in tqdm(reads_lst):
	sub_strings = [ref_read[i:i+k] for i in range(len(ref_read)-k)]
	for sub in sub_strings:
		kmer_dict.setdefault(sub,0)
		kmer_dict[sub] += 1

ps = np.array(list(kmer_dict.values()))
ps = ps/ps.sum()


print("parallelizing randRead generation")
pool      = mp.Pool(processes=num_jobs)
arg_tuple =  itertools.product(range(numReads), [dataset], [numHashes])
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=numReads):
	pass
