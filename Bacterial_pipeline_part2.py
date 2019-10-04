#!/usr/bin/env python
# coding: utf-8

# In[3]:


import sys
import os
import operator
import numpy as np
import pickle
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO
import json
import subprocess


# In[1]:


base_path = "/data/MAB_alignment/"

bact_name = "positive_examples"

if len(sys.argv) >= 3:
    base_path = sys.argv[1]
    bact_name = sys.argv[2]


st_point = 0 

if len(sys.argv) >= 4:
    st_point = int(sys.argv[3])



cwd_path = base_path+bact_name+"_filtered/"

print(cwd_path)


# In[4]:





# In[8]:


if st_point <= 1:
	subprocess.call("rm -f *.db",shell=True,cwd=cwd_path)
	fasta2DB_cmd = "fasta2DB "+bact_name+" reads.fasta"
	print (fasta2DB_cmd)
	subprocess.check_output(fasta2DB_cmd.split(),cwd=cwd_path)


# In[21]:


if st_point <= 2:
    subprocess.call("rm -f *.las",shell=True,cwd=cwd_path)
    daligner_cmd = "HPC.daligner "+bact_name
    daligner_shell_cmd = "csh -v daligner_cmd.sh"
    print (daligner_cmd)
    p = subprocess.call(daligner_cmd.split(),stdout=open(cwd_path+'daligner_cmd.sh','w') , cwd=cwd_path)
    p2 = subprocess.check_output(daligner_shell_cmd.split(), cwd=cwd_path)


# In[25]:


if st_point <= 3:
    aln_file = cwd_path+"alignments.txt"
    LAshow_cmd = "LAshow "+bact_name+" "+bact_name
    print(LAshow_cmd)
    subprocess.call(LAshow_cmd.split(),stdout=open(aln_file,'w'),cwd=cwd_path)


# In[22]:


if st_point <= 4:
    fasta_name = cwd_path+"reads.fasta"
    rd_lgths = {}
    fasta_sequences = SeqIO.parse(open(fasta_name),'fasta')

    for ind, fasta in enumerate(fasta_sequences):
        rd_lgths[ind] = len(str(fasta.seq))


# In[27]:


if st_point <= 5:
    output_txt = cwd_path+"daligner_ground_truth.txt"
    read_no = 1
    lst_alignments = []
    aln_size = {}
    reads_to_keep = set()

    with open(output_txt,'w') as f:
        pass

    with open(aln_file, 'r') as f:
        for ind, line in enumerate(f):
            if ind < 2:
                continue
            x = line.split()
            y = line.replace('<', '[').split("[")
    #         print(y)
            try:
                z = y[1].strip().split("x")[0].strip().strip(">").strip("]").split("..")
            except:
                print (line, x, y)
                raise
            start_pos = int(z[0].strip().replace(',', '')) 
            end_pos = int(z[1].strip().replace(',', ''))
    #         print len(x)


            if int(x [0].replace(',', '')) == read_no:
                if x [2] == 'n':
                    read_b_n0 = int(x [1].replace(',', ''))
                    raw_aln_size = (end_pos-start_pos+0.0)
                    aln = raw_aln_size
    #                 print read_b_n0, start_pos, end_pos, aln
                    if read_b_n0 in aln_size:
                        aln_size[read_b_n0] += aln
                    else:
                        aln_size[read_b_n0] = aln
    #                 print x
            else:
    #             sorted_aln = sorted(aln_size.items(), key=operator.itemgetter(1),reverse=True)
    #             if rd_lgths[read_no-1] > 7000 and len(sorted_aln) > 5 and sorted_aln[4][1] > 0.8:

    #                 reads_to_keep = reads_to_keep.union(set([read_no]))
    #                 reads_to_keep = reads_to_keep.union(set([sorted_aln[i][0] for i in range(5)]))
                with open(output_txt,'a') as f:
                    for rd_nm in aln_size:
                        ln_write = str(read_no)+"\t"+str(rd_nm)+"\t"+str(aln_size[rd_nm])+"\t"+str(rd_lgths[read_no-1])+"\t"+str(rd_lgths[rd_nm-1])+"\n"
                        f.write(ln_write)

                read_no  = int(x [0].replace(',', ''))
                aln_size = {}
                if x [2] == 'n':
                    read_b_n0 = int(x [1].replace(',', ''))
                    raw_aln_size = (end_pos-start_pos+0.0)
                    aln = raw_aln_size
    #                 print read_b_n0, start_pos, end_pos, aln
                    if read_b_n0 in aln_size:
                        aln_size[read_b_n0] += aln
                    else:
                        aln_size[read_b_n0] = aln


# In[ ]:




