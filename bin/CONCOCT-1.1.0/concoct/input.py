import re
import logging

import numpy as np
import pandas as p

from itertools import product, tee
from collections import Counter, OrderedDict

from Bio import SeqIO

def load_data(args):
    composition, contig_lengths = load_composition(
        args.composition_file, 
        args.kmer_length, 
        args.length_threshold
        )

    if args.coverage_file:
        cov, cov_range = load_coverage(
            args.coverage_file,
            contig_lengths,
            args.no_cov_normalization,
            add_total_coverage = (not args.no_total_coverage),
            read_length = args.read_length
            )
    else:
        cov, cov_range = None, None

    return composition, cov, cov_range


def _calculate_composition(comp_file, length_threshold, kmer_len):
    #Generate kmer dictionary
    feature_mapping, nr_features = generate_feature_mapping(kmer_len)

    # Store composition vectors in a dictionary before creating dataframe
    composition_d = OrderedDict()
    contig_lengths = OrderedDict()
    for seq in SeqIO.parse(comp_file,"fasta"):
        seq_len = len(seq)
        if seq_len<= length_threshold:
            continue
        contig_lengths[seq.id] = seq_len
        # Create a list containing all kmers, translated to integers
        kmers = [
                feature_mapping[kmer_tuple]
                for kmer_tuple 
                in window(str(seq.seq).upper(), kmer_len)
                if kmer_tuple in feature_mapping
                ]
        # numpy.bincount returns an array of size = max + 1
        # so we add the max value and remove it afterwards
        # numpy.bincount was found to be much more efficient than
        # counting manually or using collections.Counter
        kmers.append(nr_features - 1)
        composition_v = np.bincount(np.array(kmers))
        composition_v[-1] -= 1
        # Adding pseudo counts before storing in dict
        composition_d[seq.id] = composition_v + np.ones(nr_features)
    composition = p.DataFrame.from_dict(composition_d, orient='index', dtype=float)
    contig_lengths = p.Series(contig_lengths, dtype=float)
    return composition, contig_lengths

def load_composition(comp_file, kmer_len, threshold):
    #Composition
    composition, contig_lengths = _calculate_composition(
            comp_file, 
            threshold,
            kmer_len
            )

    #Normalize kmer frequencies to remove effect of contig length
    #log(p_ij) = log[(X_ij +1) / rowSum(X_ij+1)]
    composition = np.log(composition.divide(composition.sum(axis=1),axis=0))
    
    logging.info('Successfully loaded composition data.')
    return composition, contig_lengths

def load_coverage(cov_file, contig_lengths, no_cov_normalization, add_total_coverage=False, read_length=100):
    #Coverage import, file has header and contig ids as index
    cov = p.read_table(cov_file, header=0, index_col=0)

    cov = cov[cov.index.isin(contig_lengths.index)]

    # cov_range variable left here for historical reasons. Can be removed entirely
    cov_range = (cov.columns[0],cov.columns[-1])

    # Adding pseudo count
    cov.loc[:,cov_range[0]:cov_range[1]] = cov.loc[:,cov_range[0]:cov_range[1]].add(
            (read_length/contig_lengths),
            axis='index')

    if not no_cov_normalization:
        #Normalize per sample first
        cov.loc[:,cov_range[0]:cov_range[1]] = \
            _normalize_per_sample(cov.loc[:,cov_range[0]:cov_range[1]])

    temp_cov_range = None
    # Total coverage should be calculated after per sample normalization
    if add_total_coverage:
        cov['total_coverage'] = cov.loc[:,cov_range[0]:cov_range[1]].sum(axis=1)
        temp_cov_range = (cov_range[0],'total_coverage')
    
    if not no_cov_normalization:
        # Normalize contigs next
        cov.loc[:,cov_range[0]:cov_range[1]] = \
            _normalize_per_contig(cov.loc[:,cov_range[0]:cov_range[1]])

    if temp_cov_range:
        cov_range = temp_cov_range

    # Log transform
    cov.loc[:,cov_range[0]:cov_range[1]] = np.log(
        cov.loc[:,cov_range[0]:cov_range[1]])

    logging.info('Successfully loaded coverage data.')
    return cov, cov_range
    
def _normalize_per_sample(arr):
    """ Divides respective column of arr with its sum. """
    return arr.divide(arr.sum(axis=0),axis=1)

def _normalize_per_contig(arr):
    """ Divides respective row of arr with its sum. """
    return arr.divide(arr.sum(axis=1),axis=0)
    

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)
