import os

cmd_0 = "python makeFolders.py"
print(cmd_0)
os.system(cmd_0)

datasetLst = [line.rstrip('\n') for line in open('NCTC_ds.txt')]
for i in range(len(datasetLst)):
    if len(datasetLst[i])==0:
        break
    bact = datasetLst[i]

    cmd_1 = "python generateMinHashes_pipelined.py --dataset "+bact+"_filtered --num_jobs 50"
    print(cmd_1)
    os.system(cmd_1)
    
    cmd_12 = "python generateRandReads.py --dataset "+bact+"_filtered --num_jobs 50 --num_reads 5"
    print(cmd_12)
    os.system(cmd_12)

    cmd_21 = "python computeJSim.py --dataset "+bact+"_filtered --num_jobs 50"
    print(cmd_21)
    os.system(cmd_21)
    
    cmd_2 = "python runSVDexperiments_pipelined.py --dataset "+bact+"_filtered --num_jobs 50 --num_rand 5"
    print(cmd_2)
    os.system(cmd_2)
    
cmd_3 = "python generateAUCs_pipelined.py --datasets NCTC_ds.txt --num_jobs 50 --ground_truth daligner --theta .3"
print(cmd_3)
os.system(cmd_3)

cmd_4 = "python generateColProb_pipelined.py --datasets NCTC_ds.txt --num_jobs 50"
print(cmd_4)
os.system(cmd_4)
