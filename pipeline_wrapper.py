import os
import argparse

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--datasets",  help="Text file with folder to dataset on each line", type = str, default = "NCTC_ds.txt")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=50 )
ap.add_argument("--ground_truth", help="Which aligner ground truth to use (e.g. minimap2,daligner) ", type = str, default = "daligner")
ap.add_argument("--theta", help="Alignment threshold to detect, decimal between 0 and 1", type = float, default=0.3)
ap.add_argument("--num_rand", help="Number of random calibration reads", type = int, default=5)

args = ap.parse_args()

num_jobs   = args.num_jobs
datasets   = args.datasets
ground_truth = args.ground_truth
theta = args.theta
numRandReads = args.num_rand


cmd_0 = "python makeFolders.py --datasets {}".format(datasets)
print(cmd_0)
# os.system(cmd_0)

datasetLst = [line.rstrip('\n') for line in open(datasets)]
for i in range(len(datasetLst)):
    if len(datasetLst[i])==0:
        break
    bact = datasetLst[i]

    cmd_1 = "python generateMinHashes_pipelined.py --dataset {}_filtered --num_jobs {}".format(
    bact, num_jobs)
    print(cmd_1)
    os.system(cmd_1)

    cmd_12 = "python generateRandReads.py  --dataset {}_filtered --num_jobs {}  --num_reads {}".format(
    	bact, num_jobs, numRandReads)
    print(cmd_12)
    os.system(cmd_12)

    cmd_21 = "python computeJSim.py --dataset {}_filtered --num_jobs {}".format(
    bact, num_jobs)
    print(cmd_21)
    os.system(cmd_21)

    cmd_2 = "python runSVDexperiments_pipelined.py --dataset {}_filtered --num_jobs {}  --num_reads {}".format(
    	bact, num_jobs, numRandReads)
    print(cmd_2)
    os.system(cmd_2)
    # break

cmd_3 = "python generateAUCs_pipelined.py --datasets {} --num_jobs {} --ground_truth {} --theta {}".format(datasets,
num_jobs, ground_truth, theta)
print(cmd_3)
os.system(cmd_3)

cmd_4 = "python generateColProb_pipelined.py --datasets {} --num_jobs {}".format(
datasets, num_jobs)
print(cmd_4)
os.system(cmd_4)
