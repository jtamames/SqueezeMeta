import gzip

class SeqReader:

    # Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
    def readfq(self, fp): # this is a generator function
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break

    def read_nucleotide_sequences(self, nucleotide_file):
        nucleotide_sequences = {}
        for name, seq, _ in self.readfq(open(nucleotide_file)):
            nucleotide_sequences[name] = seq
        return nucleotide_sequences


    # def check_for_proper_nucleotide_seq(self, seq_file, req_perc=0.9, max_seqs_to_read=10):
    #
    #
    #     """Check if a file contains sequences in nucleotide space.
    #     The check is performed by looking for the characters in
    #     {a,c,g,t,n,.,-} and confirming that these comprise the
    #     majority of a sequences. A set number of sequences are
    #     read and the file assumed to be not be in nucleotide space
    #     if none of these sequences are comprised primarily of the
    #     defined nucleotide set.
    #     Parameters
    #     ----------
    #     seq_file : str
    #         Name of fasta/q file to read.
    #     req_perc : float
    #         Percentage of bases in {a,c,g,t,n,.,-} before
    #         declaring the sequences as being in nucleotide
    #         space.
    #     max_seqs_to_read : int
    #         Maximum sequences to read before declaring
    #         sequence file to not be in nucleotide space.
    #     Returns
    #     -------
    #     boolean
    #         True is sequences are in nucleotide space, or file
    #         contains no sequences.
    #     """
    #
    #     nucleotide_bases = {'a', 'c', 'g', 't'}
    #     insertion_bases = {'-', '.'}
    #
    #     seqs = self.read_nucleotide_sequences(seq_file)
    #     if len(seqs) == 0:
    #         return True
    #
    #     seq_count = 0
    #     for _seq_id, seq in seqs.items():
    #         seq = seq.lower()
    #
    #         nt_bases = 0
    #         for c in (nucleotide_bases | {'n'} | insertion_bases):
    #             nt_bases += seq.count(c)
    #
    #         if float(nt_bases) / len(seq) >= req_perc:
    #             return True
    #
    #         seq_count += 1
    #         if seq_count == max_seqs_to_read:
    #             break
    #
    #     return False

    def write_fasta(self, seq, outputFile):
        '''write sequences to FASTA file'''
        if outputFile.endswith('.gz'):
            fout = gzip.open(outputFile, 'wb')
        else:
            fout = open(outputFile, 'w')

        for seqId, seq in seq.items():
            fout.write('>' + seqId + '\n')
            fout.write(seq + '\n')
        fout.close()