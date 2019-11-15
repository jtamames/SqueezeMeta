###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from collections import defaultdict

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.common import remove_extension

"""
Functions for verify, exploring, modifying and
calculating statistics on one or more genomes.
"""


def gc_count(seqs):
    """Determine number of G or C bases in sequences.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        Number of G or C bases in sequences.
    """

    C = 0
    G = 0
    for seq in list(seqs.values()):
        _a, c, g, _t = seq_tk.count_nt(seq)

        C += c
        G += g

    return G + C


def gc(seqs):
    """Calculate GC of sequences.

    GC is calculated as (G+C)/(A+C+G+T), where
    each of these terms represents the number
    of nucleotides within the sequence. Ambiguous
    and degenerate bases are ignored. Uracil (U)
    is treated as a thymine (T).

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    float
        GC content of sequences.
    """

    A = 0
    C = 0
    G = 0
    T = 0
    for seq in list(seqs.values()):
        a, c, g, t = seq_tk.count_nt(seq)

        A += a
        C += c
        G += g
        T += t

    return float(G + C) / (A + C + G + T)


def ambiguous_nucleotides(seqs):
    """Count ambiguous or degenerate nucleotides in sequences.

    Any base that is not a A, C, G, or T/U is considered
    to be ambiguous or degenerate.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        Number of ambiguous and degenerate bases.
    """

    num_ambiguous = 0
    for seq in list(seqs.values()):
        num_ambiguous += seq_tk.ambiguous_nucleotides(seq)

    return num_ambiguous


def unique(genome_files):
    """Check if sequences are assigned to multiple bins.

    Parameters
    ----------
    genome_files : iterable
        Path to genome fasta files.

    Returns
    -------
    dict : d[genome_id][genome_id] -> [shared sequences]
        List of any sequences within a genome observed multiple times.
    """

    # read sequence IDs from all genomes,
    # while checking for duplicate sequences within a genomes
    duplicates = defaultdict(lambda: defaultdict(list))

    genome_seqs = {}
    for f in genome_files:
        genome_id = remove_extension(f)

        seq_ids = set()
        for seq_id, _seq in seq_io.read_seq(f):
            if seq_id in seq_ids:
                duplicates[genome_id][genome_id].append(seq_id)

            seq_ids.add(seq_id)

        genome_seqs[genome_id] = seq_ids

    # check for sequences assigned to multiple bins
    genome_ids = list(genome_seqs.keys())
    for i in range(0, len(genome_ids)):
        seq_idsI = genome_seqs[genome_ids[i]]

        for j in range(i + 1, len(genome_ids)):
            seq_idsJ = genome_seqs[genome_ids[j]]

            seq_intersection = seq_idsI.intersection(seq_idsJ)

            if len(seq_intersection) > 0:
                duplicates[genome_ids[i]][genome_ids[j]] = seq_intersection
                duplicates[genome_ids[j]][genome_ids[i]] = seq_intersection

    return duplicates


def modify(input_file, scaffold_file, seqs_to_add, seqs_to_remove, output_file):
    """Add or remove scaffolds from a fasta file.

    Parameters
    ----------
    input_file : str
        Fasta file to modify.
    scaffold_file : str
        Fasta file containing scaffolds to add.
    seqs_to_add: iterable
        Unique ids of scaffolds to add.
    seqs_to_remove : iterable
        Unique ids of scaffolds to remove.
    output_file : str
        Desired name of modified fasta file.

    Returns
    -------
    iterable, iterable
        Unique ids of sequences that could not be added,
        unique ids of sequences that could not be removed.
    """

    seqs = seq_io.read(input_file)

    # add sequences to bin
    failed_to_add = set()
    if seqs_to_add:
        failed_to_add = set(seqs_to_add)
        if seqs_to_add != None:
            for seq_id, seq in seq_io.read_seq(scaffold_file):
                if seq_id in seqs_to_add:
                    failed_to_add.remove(seq_id)
                    seqs[seq_id] = seq

    # remove sequences from bin
    failed_to_remove = set()
    if seqs_to_remove:
        failed_to_remove = set(seqs_to_remove)
        if seqs_to_remove != None:
            for seq_id in seqs_to_remove:
                if seq_id in seqs:
                    failed_to_remove.remove(seq_id)
                    seqs.pop(seq_id)

    # save modified bin
    seq_io.write_fasta(seqs, output_file)

    return failed_to_add, failed_to_remove
