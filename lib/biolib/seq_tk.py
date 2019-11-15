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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'


#import string # Commented by FPS NOV 15 2019, we now use str.maketrans instead
from collections import Counter


"""Sequence manipulation and statistics."""

_complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # FPS NOV 15 2019


def count_nt(seq):
    """Count occurrences of each nucleotide in a sequence.

    Only the bases A, C, G, and T(U) are counted. Ambiguous
    bases are ignored.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    list
        Number of A, C, G, and T(U) in sequence.
    """

    s = seq.upper()
    a = s.count('A')
    c = s.count('C')
    g = s.count('G')
    t = s.count('T') + s.count('U')

    return a, c, g, t


def gc(seq):
    """Calculate GC content of a sequence.

    GC is calculated as (G+C)/(A+C+G+T), where
    each of these terms represents the number
    of nucleotides within the sequence. Ambiguous
    and degenerate bases are ignored. Uracil (U)
    is treated as a thymine (T).

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    float
        GC content of sequence.
    """

    a, c, g, t = count_nt(seq)
    return float(g + c) / (a + c + g + t)


def ambiguous_nucleotides(seq):
    """Count ambiguous or degenerate nucleotides in a sequence.

    Any base that is not a A, C, G, or T/U is considered
    to be ambiguous or degenerate.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    int
        Number of ambiguous and degenerate bases.
    """

    a, c, g, t = count_nt(seq)
    return len(seq) - (a + c + g + t)


def rev_comp(seq):
    """Reverse complement a sequence."""
    return seq.translate(_complements)[::-1]


def N50(seqs):
    """Calculate N50 for a set of sequences.

     N50 is defined as the length of the longest
     sequence, L, for which 50% of the total bases
     are present in sequences of length >= L.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        N50 for the set of sequences.
    """

    seq_lens = [len(x) for x in list(seqs.values())]
    threshold = sum(seq_lens) / 2.0

    seq_lens.sort(reverse=True)

    current_sum = 0
    for seq_len in seq_lens:
        current_sum += seq_len
        if current_sum >= threshold:
            N50 = seq_len
            break

    return N50


def L50(seqs, N50):
    """Calculate L50 for a set of sequences.

     L50 is defined as the number of sequences
     that are longer than, or equal to, the
     N50 length.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.
    N50 : int
        N50 of sequences.

    Returns
    -------
    int
        L50 for the set of sequences.
    """

    L50 = sum([1 for x in list(seqs.values()) if len(x) >= N50])

    return L50


def mean_length(seqs):
    """Calculate mean length of sequences.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    float
        Mean length of sequences.
    """

    total_len = sum([len(x) for x in list(seqs.values())])

    return float(total_len) / len(seqs)


def max_length(seqs):
    """Identify longest sequence.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        Length of longest sequence.
    """

    return max([len(x) for x in list(seqs.values())])


def identify_contigs(seqs, contig_break='NNNNNNNNNN'):
    """Break scaffolds into contigs.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.
    contig_break : str
        Motif used to split scaffolds into contigs.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Contigs indexed by sequence ids.
    """

    contigs = {}
    for seq_id, seq in seqs.items():
        seq = seq.upper()
        contig_count = 0
        for contig in seq.split(contig_break):
            contig = contig.replace('N', '')
            if contig:
                contigs[seq_id + '_c' + str(contig_count)] = contig
                contig_count += 1

    return contigs


def fragment(seq, window_size, step_size):
    """Fragment sequence into fixed sized windows.

    The last fragment may not be shorter than
    the window size, but will only be generated
    if it is at least half the window size.

    Parameters
    ----------
    seq : str
        Sequence to fragment.
    window_size : int
        Size of each fragment.
    step_size : int
        Number of bases to move after each window.

    Returns
    -------
    list
        Fragments from sequences.
    """

    fragments = []
    start = 0
    for i in range(0, len(seq), step_size):
        end = i + window_size
        if end < len(seq):
            fragments.append(seq[start:end])
            start = end

    # get last fragment if it is at least half
    # the specified window size
    if len(seq) - start >= 0.5 * window_size:
        fragments.append(seq[start:])

    return fragments


def trim_seqs(seqs, min_per_taxa, consensus, min_per_bp):
        """Trim multiple sequence alignment.

        Adapted from the biolib package.

        Parameters
        ----------
        seqs : d[seq_id] -> sequence
            Aligned sequences.
        min_per_taxa : float
            Minimum percentage of taxa required to retain a column [0,1].
        min_per_bp : float
            Minimum percentage of base pairs required to keep trimmed sequence [0,1].
        Returns
        -------
        dict : d[seq_id] -> sequence
            Dictionary of trimmed sequences.
        dict : d[seq_id] -> sequence
            Dictionary of pruned sequences.
        int 
            Number of columns filtered by minimum percentage of taxa.
        int 
            Number of columns filtered by consensus
        """

        alignment_length = len(list(seqs.values())[0])

        # count number of taxa represented in each column
        column_count = [0] * alignment_length
        column_chars = [list() for _ in range(alignment_length)]
        for seq in list(seqs.values()):
            for i, ch in enumerate(seq):
                if ch != '.' and ch != '-':
                    column_count[i] += 1
                    column_chars[i].append(ch)

        mask = [False] * alignment_length
        count_min_taxa_filtered = 0
        count_consensus_filtered = 0
        for i, count in enumerate(column_count):
            if count >= min_per_taxa * len(seqs):
                c = Counter(column_chars[i])
                if len(c.most_common(1)) == 0:
                    ratio = 0
                else:
                    _letter, count = c.most_common(1)[0]
                    ratio = float(count) / column_count[i]
                    
                if ratio >= consensus:
                    mask[i] = True
                else:
                    count_consensus_filtered += 1
            else:
                count_min_taxa_filtered += 1

        # trim columns
        output_seqs = {}
        pruned_seqs = {}
        for seq_id, seq in seqs.items():
            masked_seq = ''.join([seq[i] for i in range(0, len(mask)) if mask[i]])

            valid_bases = len(masked_seq) - masked_seq.count('.') - masked_seq.count('-')
            if valid_bases < len(masked_seq) * min_per_bp:
                pruned_seqs[seq_id] = masked_seq
                continue

            output_seqs[seq_id] = masked_seq

        return output_seqs, pruned_seqs, count_min_taxa_filtered, count_consensus_filtered


def aai(seq1, seq2):
    """Calculate AAI between sequences.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : float
        Second sequence.

    Returns
    -------
    float
        AAI between sequences.
    """

    assert len(seq1) == len(seq2)

    mismatches = 0
    matches = 0
    for c1, c2 in zip(seq1.upper(), seq2.upper()):
        if c1 != c2:
            mismatches += 1
        elif c1 == '-' and c2 == '-':
            pass
        else:
            matches += 1

        return matches * 100.0 / (matches + mismatches)
