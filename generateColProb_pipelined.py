from Bio.SeqIO.FastaIO import SimpleFastaParser
import mmh3
import itertools
import multiprocessing as mp
import numpy as np
import os, pickle
import argparse
from tqdm import tqdm, trange


def runSim(args):

	## avoid one processes starting multiple threads
	os.environ["MKL_NUM_THREADS"] = "1" 
	os.environ["NUMEXPR_NUM_THREADS"] = "1" 
	os.environ["OMP_NUM_THREADS"] = "1" 

	dataset = args
	if dataset[-1] != '/':
		dataset += '/'

	## load dna reads into reads_lst 
	reads_lst = []
	fastaFile = dataset + "/reads.fasta"
	with open(fastaFile) as handle:
		for values in SimpleFastaParser(handle):
			reads_lst.append(values[1])
	n = len(reads_lst)

	symLength=7 # k in k-mer

	def generateSymSets(reads_lst,symLength):
		symSets = {}
		for i,read in enumerate(reads_lst):
			lst = set()
			for j in range(len(read)-symLength):
				lst.add(read[j:j+symLength])
			symSets[i] = lst
		return symSets

	symSets = generateSymSets(reads_lst,symLength)


	## Calculate collision probability of kmer dist
	kmer_dict = {}
	k = 7
	for ref_read in reads_lst:
		sub_strings = [ref_read[i:i+k] for i in range(len(ref_read)-k)]
		for sub in sub_strings:
			kmer_dict.setdefault(sub,0)
			kmer_dict[sub] += 1

	normalization = np.array(list(kmer_dict.values())).sum()
	for k in kmer_dict:
		kmer_dict[k] = kmer_dict[k]/normalization

	mapping = ['a','t','c','g']
	def intToDNA(x):
		ret = ''
		for _ in range(7):
			ret += mapping[x%4]
			x//=4
		return ret

	numHashes = 1000

	hashPerm = np.zeros((numHashes,4**7))
	for hashNum in trange(numHashes):
		tmp = [mmh3.hash(intToDNA(i),hashNum,signed=False) for i in range(4**7)]
		hashPerm[hashNum] = np.argsort(np.array(tmp)).astype(int)

	colProbStoreArrMean = []
	colProbStoreArrMedian = []
	colProbStoreArrAll = []
	argminKmerProb = []

	for ref_read_idx in trange(n):

		colProbArr = []
		ref_read_sym = symSets[ref_read_idx]
		L=int(np.mean([len(x) for x in reads_lst]))

		for hashNum in range(1000):
			idx = hashPerm[hashNum,:].astype(int)
			argminKmerRank = 0
			for i,dnaVal in enumerate(idx):
				if intToDNA(int(dnaVal)) in ref_read_sym:
					argminKmerRank = i
					break

			## numerically stable approximation to collision probability
			pre = np.sum([kmer_dict[intToDNA(x)] for x in idx[:argminKmerRank]])
			prob = np.exp(-L*pre) - np.exp(-L*(pre+kmer_dict[intToDNA(idx[argminKmerRank])]))
	
			colProbArr.append(prob)
			argminKmerProb.append(kmer_dict[intToDNA(idx[argminKmerRank])])
		
		colProbStoreArrAll.extend(colProbArr)
		colProbStoreArrMean.append(np.mean(colProbArr))
		colProbStoreArrMedian.append(np.median(colProbArr))

	pickle.dump([colProbStoreArrAll,colProbStoreArrMean,argminKmerProb,
		"colProbStoreArrAll, colProbStoreArrMean, argminKmerProb"],
		open("AUCs/{}_colProb.pkl".format(dataset[:-1]), "wb" ) )


### parallelizing
ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--datasets",  help="Text file with folder to dataset on each line", type = str, default = "NCTC_ds.txt")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )

args = ap.parse_args()
num_jobs   = args.num_jobs
datasets   = args.datasets      

if not os.path.isdir("AUCs"):
	os.mkdir("AUCs")

datasetLst = [line.rstrip('\n') for line in open(datasets)]
for i in range(len(datasetLst)):
	if len(datasetLst[i])==0:
		datasetLst = datasetLst[:i]
		break

datasetLst = [x+"_filtered" for x in datasetLst]
print("running on ",datasetLst)
num_jobs = min(num_jobs,len(datasetLst))
print("Parallelizing {} collision probability computations over {} workers".format(len(datasetLst),num_jobs))
pool      = mp.Pool(processes=num_jobs)
arg_tuple = datasetLst
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=len(datasetLst)):
	pass
