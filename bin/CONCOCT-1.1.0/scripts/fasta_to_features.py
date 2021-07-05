#!/usr/bin/env python
import numpy as np
import pandas as pd
from itertools import product
from Bio import SeqIO

# optimized sliding window function from
# http://stackoverflow.com/a/7636587
from itertools import tee

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        kmer = ''.join(kmer)
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = ''.join([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash

def generate_features_from_fasta(fasta_file,nr_datapoints,kmer_len,outfile):
    kmer_dict = generate_feature_mapping(kmer_len)
    seqs = SeqIO.parse(fasta_file,"fasta")
    #Initialize feature vectors. NxD where N is number of datapoints, D is number of dimentions
    contigs = np.zeros((nr_datapoints,max(kmer_dict.values())+1))
    contigs_id = []
    for i,seq in enumerate(seqs):
        contigs_id.append(seq.id)
        for kmer_tuple in window(str(seq.seq).upper(),kmer_len):
            contigs[i,kmer_dict["".join(kmer_tuple)]] += 1
    df = pd.DataFrame(contigs,index=contigs_id)
    df.to_csv(outfile)
#    print contigs_id
#    contigs_id = np.array(contigs_id).reshape((-1,1))
#    print contigs_id.shape
#    print contigs.shape

#    np.savetxt(outfile,np.hstack((contigs_id,contigs)),delimiter=",")
    
if __name__=="__main__":
    import sys
    fasta_file = sys.argv[1]
    nr_datapoints = int(sys.argv[2])
    kmer_len = int(sys.argv[3])
    outfile = sys.argv[4]
    generate_features_from_fasta(fasta_file,nr_datapoints,kmer_len,outfile)
