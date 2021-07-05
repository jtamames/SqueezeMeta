#!/usr/bin/env python
from nose.tools import assert_equal, assert_true
import numpy as np
import pandas as p
import os
from Bio import SeqIO
from concoct.input import _normalize_per_sample, _normalize_per_contig, generate_feature_mapping, load_composition, _calculate_composition

class TestInput(object):
    def setUp(self):
        self.C = p.DataFrame(np.array([[0., 0.7], [5.5, .7]]))

    def test_normalize_per_contig(self):
        C_norm = _normalize_per_contig(self.C)
        
        C_correct = p.DataFrame(np.array([[0., 1.],[5.5/6.2, 0.7/6.2]]))
        assert_true(np.linalg.norm(C_norm-C_correct) < 0.0001)

    def test_normalize_per_samples(self):
        C_norm = _normalize_per_sample(self.C)

        C_correct = p.DataFrame(np.array([[0., 0.5],[1,0.5]]))
        assert_true(np.linalg.norm(C_norm-C_correct) < 0.0001)

    def test_generate_feature_mapping(self):
        feature_mapping, counter = generate_feature_mapping(2)
        assert_equal(counter, 10)
        assert_equal(len(list(feature_mapping.keys())), 16)
        assert_true(('A', 'A') in feature_mapping)

    def test_load_composition(self):
        # Get the directory path of this test file
        f = os.path.dirname(os.path.abspath(__file__))
        # calculate the lengths of the contigs
        seqs = SeqIO.parse("{0}/test_data/composition_some_shortened.fa".format(f),"fasta")
        ids = []
        lengths = []
        for s in seqs:
            if len(s) <= 1000:
                continue
            ids.append(s.id)
            lengths.append(len(s))
        c_len = p.Series(lengths,index=ids,dtype=float)
        # Use load_composition to calculate contig lengths
        composition, contig_lengths = load_composition("{0}/test_data/composition_some_shortened.fa".format(f),4,1000)
        assert_equal(len(c_len), len(contig_lengths))
        # All equal
        for ix in ids:
            assert_equal(c_len.loc[ix], contig_lengths.loc[ix])
        

    def test__calculate_composition(self):
        d = os.path.dirname(os.path.abspath(__file__))
        f = "{0}/test_data/composition_some_shortened.fa".format(d)
        seqs = SeqIO.parse(f, "fasta")   

        feature_mapping, counter = generate_feature_mapping(4)
       
        seq_strings = {}
        for i, s in enumerate(seqs):
            seq_strings[s.id] = str(s.seq).upper()

        composition, contig_lengths = _calculate_composition(f, 0, 4)  
        
        # Make sure the count is correct for one specific kmer 
        kmer_s = ('A', 'C', 'G', 'T')

        for seq_id, s in seq_strings.items():
            c = count_substrings(s, "".join(kmer_s))
            assert_equal(composition.loc[seq_id, feature_mapping[kmer_s]], c+1)

        # Check that non palindromic kmers works as well:
        kmer_s = ('A', 'G', 'G', 'G')
        reverse_kmer_s = ('C', 'C', 'C', 'T')
        for seq_id, s in seq_strings.items():
            c_1 = count_substrings(s, "".join(kmer_s))
            c_2 = count_substrings(s, "".join(reverse_kmer_s))
            assert_equal(composition.loc[seq_id, feature_mapping[kmer_s]], c_1 + c_2 + 1)
        

def count_substrings(s, subs):
    # stolen from http://stackoverflow.com/a/19848382
    # modified to count overlapping substrings as well
    start = numBobs = 0 
    while start >= 0:
        pos = s.find(subs, start)
        if pos < 0:
            break
        numBobs += 1
        start = pos + 1
    return numBobs
