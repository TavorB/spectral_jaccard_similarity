from Bio.SeqIO.FastaIO import SimpleFastaParser
import mmh3
import itertools
import multiprocessing as mp
import numpy as np
import os, pickle
import argparse
from tqdm import tqdm, trange
from sklearn.metrics import roc_curve, auc



def runSim(args):

	## avoid one processes starting multiple threads
	os.environ["MKL_NUM_THREADS"] = "1" 
	os.environ["NUMEXPR_NUM_THREADS"] = "1" 
	os.environ["OMP_NUM_THREADS"] = "1" 

	dataset = args[0]
	if dataset[-1] != '/':
		dataset += '/'
	ground_truth,theta = args[1]

	## load dna reads into reads_lst 
	reads_lst = []
	fastaFile = dataset + "/reads.fasta"
	with open(fastaFile) as handle:
		for values in SimpleFastaParser(handle):
			reads_lst.append(values[1])
	n = len(reads_lst)

	## load precomputed Jaccard Similarities
	JSims = np.loadtxt(dataset+"/minHashes/JSims.txt")


	## load alignments
	gt_file = "{}/{}_ground_truth.txt".format(dataset,ground_truth)

	with open(gt_file) as f:
		lines = [[float(x) for x in line.rstrip('\n').split('\t')] for line in f]

	refDict = {}
	for i in range(n):
		refDict[i] = {}
	for line in lines:
		refDict[int(line[0])-1][int(line[1])-1] = line[2]/(line[3]+line[4]-line[2])
		refDict[int(line[1])-1][int(line[0])-1] = line[2]/(line[3]+line[4]-line[2])


	## convert each read into it's k-mers
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

        
	## load precomputed minHashes
	minHashArr = np.zeros((n,1000))
	for i in range(n):
		minHashArr[i] = np.load(dataset+"minHashes/minHashes_{}.txt".format(i), allow_pickle=True)

	## load precomputed iid sequences and their minhashes
	numRandReads = 5
	randMinHashArr = np.zeros((numRandReads,1000))
	for i in range(numRandReads):
		randMinHashArr[i] = np.load(dataset+"randReads/randMinHashes_{}.txt".format(i), allow_pickle=True)

	minHashArrExtended = np.vstack((minHashArr[:,:1000],randMinHashArr))


	## checking to see precomputed minHashes works
	i = np.random.randint(0,n,size = 100)
	j = np.random.randint(0,1000,size = 100)
	lst = []
	for iter_round in range(100):
		
		iterLst = list(symSets[i[iter_round]])
		
		lst.append(min([mmh3.hash(sym, j[iter_round], signed=False) for sym in iterLst]))
		
	assert(np.alltrue(minHashArrExtended[i,j] == lst))


	## checking to see precomputed JSim works
	i = np.random.randint(0,n,size = 100)
	j = np.random.randint(0,n,size = 100)
	lst = []
	for iter_round in range(100):
		i1 = i[iter_round]
		j1 = j[iter_round]

		lst.append(JSims[i1,j1]==1.0*len(symSets[i1].intersection(symSets[j1]))/(len(symSets[i1].union(symSets[j1]))))
		
	assert(np.alltrue(lst) and np.allclose(JSims,JSims.T))


	## Testing SVD, JSimEmp, JSim Exact, reference vs all
	storageArrGround = []
	storageArrpHatSVD = []
	storageArrJsimExact = []
	storageArrJsimEmp = []
	storageArrNumOnesCol = []
	storageArrQjs = []
	storageArrwSJS = []

	h = 1000

	for ref_read in trange(n):

		groundTruthLocs = np.array(list(refDict[ref_read].keys()))
		if len(groundTruthLocs)==0: ## read has no alignments in dataset
			continue
		refReadMatches = refDict[ref_read]
		groundTruthVals = [refReadMatches[i] for i in groundTruthLocs]

		lst = set(list(groundTruthLocs))
		rangeN = set(range(n))
		rangeN.discard(ref_read)
		
		toAppend = list(rangeN-lst)
		groundTruthLocs = np.hstack((groundTruthLocs,np.array(toAppend)))
		groundTruthVals += [0]*len(toAppend)

		empiricalMatrix = (minHashArrExtended == minHashArrExtended[ref_read])
		empiricalMatrix = np.delete(empiricalMatrix,ref_read,axis=0)[:,:h]

		updatedGroundTruthLocs = groundTruthLocs - 1*(groundTruthLocs>=ref_read)
		updatedGroundTruthLocs = updatedGroundTruthLocs.astype(int)

		jSimEmpirical = np.mean(empiricalMatrix,axis=1)
		jSimExact = np.delete(JSims[ref_read],ref_read)

		## here we can modify what normalization is used without having to rerun SVDs
		# u = np.loadtxt(dataset+"/SVD/raw_pi_refread_{}.txt".format(ref_read))
		# pHatSVD = 1-np.abs(u[:n-1])/np.abs(np.median(u[:n-1])) ## normalize median of p_i
		# pHatSVD = 1-np.abs(u[:n-1])/np.abs(np.median(u[n-1:])) ## random read normalization
		# pHatSVD = 1-np.abs(u[:n-1])/np.max(np.abs(u[n-1:])) ## naive max normalziation

		pHatSVD = np.loadtxt(dataset+"/SVD/pi_refread_{}.txt".format(ref_read))
		qSVD = np.loadtxt(dataset+"/SVD/qj_refread_{}.txt".format(ref_read))


		## for approximation
		empQ = empiricalMatrix.sum(axis=0)
		x = np.matmul(empiricalMatrix-np.ones(empiricalMatrix.shape),1-np.array(empQ/np.max(empQ)))[:n-1]
		x = np.abs(x- np.min(x)) 
		x/= np.max(x)
		storageArrwSJS.extend(x[updatedGroundTruthLocs])

		storageArrGround.extend(groundTruthVals)
		storageArrpHatSVD.extend(pHatSVD[updatedGroundTruthLocs])
		storageArrJsimEmp.extend(jSimEmpirical[:n-1][updatedGroundTruthLocs])
		storageArrJsimExact.extend(jSimExact[updatedGroundTruthLocs])
		storageArrNumOnesCol.extend(np.mean(empiricalMatrix,axis=0))
		storageArrQjs.extend(qSVD)

	fpr, tpr, _ = roc_curve(np.array(storageArrGround)>=theta,storageArrpHatSVD)
	fpr_jsim, tpr_jsim, _ = roc_curve(np.array(storageArrGround)>=theta,storageArrJsimExact)
	fpr_js_emp, tpr_js_emp, _ = roc_curve(np.array(storageArrGround)>=theta,storageArrJsimEmp)
	fpr_wsjs,tpr_wsjs,_ = roc_curve(np.array(storageArrGround)>=theta,storageArrwSJS)

	pickle.dump([auc(fpr,tpr),
		auc(fpr_jsim,tpr_jsim),
		auc(fpr_js_emp, tpr_js_emp),
		auc(fpr_wsjs, tpr_wsjs),
		storageArrNumOnesCol,
		np.corrcoef(storageArrpHatSVD,storageArrGround)[0,1],
		np.corrcoef(storageArrJsimExact,storageArrGround)[0,1],
		np.corrcoef(storageArrJsimEmp,storageArrGround)[0,1],
		storageArrQjs,
		'SJS AUC,JS AUC, JS emp AUC, wSJS AUC,numOnes per col,SJS r^2,JS r^2,JS emp r^2,storageArrQjs'],
		open("AUCs/{}_{}_{}.pkl".format(dataset[:-1],ground_truth,str(theta%1).split('.')[1]), "wb" ) )


### parallelizing
ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--datasets",  help="Text file with folder to dataset on each line", type = str, default = "NCTC_ds.txt")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--ground_truth", help="Which aligner ground truth to use (e.g. minimap2,daligner) ", type = str, default = "daligner")
ap.add_argument("--theta", help="Alignment threshold to detect, decimal between 0 and 1", type = float, default=0.3)
args = ap.parse_args()

num_jobs   = args.num_jobs
datasets   = args.datasets      
ground_truth = args.ground_truth
theta = args.theta


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
print("Parallelizing {} AUC computations over {} workers".format(len(datasetLst),num_jobs))
pool      = mp.Pool(processes=num_jobs)
arg_tuple = itertools.product(datasetLst, [[ground_truth,theta]])
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=len(datasetLst)):
	pass



