from checkm2 import sequenceClasses
import numpy as np
import pandas as pd
import os

class MetadataCalculator:

    def __init__(self, protein_file):
        self.protein_file = protein_file
        self.basename = os.path.splitext(os.path.basename(self.protein_file))[0]
        self.parsed_faa = sequenceClasses.SeqReader().read_nucleotide_sequences(self.protein_file)
#            self.parsed_genomes = [sequenceClasses.SeqReader().read_nucleotide_sequences(faa) for faa in
#                                   self.protein_files]

    def calculate_CDS(self):
        return self.basename, len(self.parsed_faa.keys())

    def calculate_amino_acid_length(self):
        return self.basename, np.array([len(v) for v in self.parsed_faa.values()]).sum()

    def calculate_amino_acid_counts(self):
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        full_aa_sequence = ''.join(value.strip() for key, value in sorted(self.parsed_faa.items()))
        aa_counts = [full_aa_sequence.count(aa) for aa in amino_acids]
        return self.basename, amino_acids, aa_counts
