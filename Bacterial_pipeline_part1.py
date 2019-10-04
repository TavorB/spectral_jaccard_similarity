#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
import argparse


# In[3]:

if (len(sys.argv) > 3):
    print("Wrong usage")

base_path = "/data/MAB_alignment/"
bact_name = "iid_genome"

if (len(sys.argv) == 3):
    base_path = sys.argv[1]
    bact_name = sys.argv[2]

num_reads_to_keep = 1000
min_read_length = 7000
num_aln = 5
aln_threshold = 0.8

org_path = base_path+bact_name+"/"
print(org_path)


# In[6]:


aln_txt_org = org_path+"alignments.txt"
all_read_fasta = org_path+bact_name+"_all.fasta"

DBshow1 = ("DBshow " + org_path + bact_name + " > " + all_read_fasta)
print(DBshow1)

LAshow1 = ("LAshow " + org_path + bact_name + " " + org_path + 
            bact_name + " > " + aln_txt_org)
print(LAshow1)

os.system(DBshow1)
os.system(LAshow1)


# In[7]:


rd_lgths = {}
fasta_sequences = SeqIO.parse(open(all_read_fasta),'fasta')

for ind, fasta in enumerate(fasta_sequences):
    rd_lgths[ind] = len(str(fasta.seq))


# In[9]:


read_no = 1
lst_alignments = []
aln_size = {}
reads_to_keep = set()




with open(aln_txt_org, 'r') as f:
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
            sorted_aln = sorted(aln_size.items(), key=operator.itemgetter(1),reverse=True)
            if (rd_lgths[read_no-1] > min_read_length 
                and len(sorted_aln) > num_aln 
                and sorted_aln[4][1] > aln_threshold):
             
                reads_to_keep = reads_to_keep.union(set([read_no]))
                reads_to_keep = reads_to_keep.union(set([sorted_aln[i][0] for i in range(5)]))
            if len(reads_to_keep) > num_reads_to_keep:
                break
#             with open(output_txt,'a') as f:
#                 for rd_nm in aln_size:
#                     ln_write = str(read_no)+"\t"+str(rd_nm)+"\t"+str(aln_size[rd_nm])+"\t"+str(rd_lgths[read_no-1])+"\t"+str(rd_lgths[rd_nm-1])+"\n"
#                     f.write(ln_write)
            
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


# In[12]:


reads_to_keep_file = org_path+"reads_filtered.txt"
np.savetxt(reads_to_keep_file,sorted(list(reads_to_keep)), fmt="%d")


# In[14]:


new_base_path = base_path+bact_name+"_filtered/"
mkdir_cmd = "mkdir -p "+new_base_path
print(mkdir_cmd)
os.system(mkdir_cmd)


# In[16]:


new_read_file = new_base_path+"reads.fasta"
DBshow_cmd2 = ("DBshow "+ org_path + bact_name + " " 
               + reads_to_keep_file + " > " + new_read_file)
print(DBshow_cmd2)
os.system(DBshow_cmd2)


# In[ ]:




