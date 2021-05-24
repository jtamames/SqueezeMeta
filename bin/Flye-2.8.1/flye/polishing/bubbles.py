#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

from __future__ import absolute_import
from __future__ import division
import logging
from bisect import bisect
from flye.six.moves import range

import multiprocessing
import signal

import flye.utils.fasta_parser as fp
import flye.config.py_cfg as cfg
from flye.polishing.alignment import shift_gaps, get_uniform_alignments
from flye.utils.sam_parser import SynchronizedSamReader
from flye.six.moves import zip


logger = logging.getLogger()


class ProfileInfo(object):
    __slots__ = ("nucl", "num_inserts", "num_deletions",
                 "num_missmatch", "coverage")

    def __init__(self):
        self.nucl = ""
        self.num_inserts = 0
        self.num_deletions = 0
        self.num_missmatch = 0
        self.coverage = 0


class Bubble(object):
    __slots__ = ("contig_id", "position", "branches", "consensus")

    def __init__(self, contig_id, position):
        self.contig_id = contig_id
        self.position = position
        self.branches = []
        self.consensus = ""


def _thread_worker(aln_reader, contigs_info, err_mode,
                   results_queue, error_queue, bubbles_file_handle,
                   bubbles_file_lock):
    """
    Will run in parallel
    """
    try:
        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

            #logger.debug("Processing {0}".format(ctg_id))
            #get top unifom alignments
            ctg_aln = get_uniform_alignments(ctg_aln, contigs_info[ctg_id].length)

            profile, aln_errors = _compute_profile(ctg_aln, err_mode,
                                                   contigs_info[ctg_id].length)
            partition, num_long_bubbles = _get_partition(profile, err_mode)
            ctg_bubbles = _get_bubble_seqs(ctg_aln, err_mode, profile, partition,
                                           contigs_info[ctg_id])
            mean_cov = sum([len(b.branches) for b in ctg_bubbles]) // (len(ctg_bubbles) + 1)
            ctg_bubbles, num_empty, num_long_branch = \
                                    _postprocess_bubbles(ctg_bubbles)
            results_queue.put((ctg_id, len(ctg_bubbles), num_long_bubbles,
                               num_empty, num_long_branch, aln_errors,
                               mean_cov))
            with bubbles_file_lock:
                _output_bubbles(ctg_bubbles, bubbles_file_handle)

            del profile
            del ctg_bubbles

    except Exception as e:
        error_queue.put(e)


def make_bubbles(alignment_path, contigs_info, contigs_path,
                 err_mode, num_proc, bubbles_out):
    """
    The main function: takes an alignment and returns bubbles
    """
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_sequence_dict(contigs_path),
                                       cfg.vals["max_read_coverage"],
                                       use_secondary=True)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    bubbles_out_lock = multiprocessing.Lock()
    bubbles_out_handle = open(bubbles_out, "w")
    for _ in range(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(aln_reader, contigs_info,
                                                     err_mode, results_queue,
                                                     error_queue, bubbles_out_handle,
                                                     bubbles_out_lock)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
            if t.exitcode == -9:
                logger.error("Looks like the system ran out of memory")
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}"
                                .format(t.exitcode))
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()
        raise

    if not error_queue.empty():
        raise error_queue.get()
    aln_reader.close()

    total_bubbles = 0
    total_long_bubbles = 0
    total_long_branches = 0
    total_empty = 0
    total_aln_errors = []
    coverage_stats = {}

    while not results_queue.empty():
        (ctg_id, num_bubbles, num_long_bubbles,
            num_empty, num_long_branch,
            aln_errors, mean_coverage) = results_queue.get()
        total_long_bubbles += num_long_bubbles
        total_long_branches += num_long_branch
        total_empty += num_empty
        total_aln_errors.extend(aln_errors)
        total_bubbles += num_bubbles
        coverage_stats[ctg_id] = mean_coverage

    mean_aln_error = sum(total_aln_errors) / (len(total_aln_errors) + 1)
    logger.debug("Generated %d bubbles", total_bubbles)
    logger.debug("Split %d long bubbles", total_long_bubbles)
    logger.debug("Skipped %d empty bubbles", total_empty)
    logger.debug("Skipped %d bubbles with long branches", total_long_branches)

    return coverage_stats, mean_aln_error


def _output_bubbles(bubbles, out_stream):
    """
    Outputs list of bubbles into file
    """
    for bubble in bubbles:
        out_stream.write(">{0} {1} {2}\n".format(bubble.contig_id,
                                        bubble.position,
                                        len(bubble.branches)))
        out_stream.write(bubble.consensus + "\n")
        for branch_id, branch in enumerate(bubble.branches):
            out_stream.write(">{0}\n".format(branch_id))
            out_stream.write(branch + "\n")

    out_stream.flush()


def _postprocess_bubbles(bubbles):
    MAX_BUBBLE = cfg.vals["max_bubble_length"]
    MAX_BRANCHES = cfg.vals["max_bubble_branches"]

    new_bubbles = []
    long_branches = 0
    empty_bubbles = 0
    for bubble in bubbles:
        if len(bubble.branches) == 0:
            #logger.debug("Empty bubble {0}".format(bubble.position))
            empty_bubbles += 1
            continue

        new_branches = []
        median_branch = (sorted(bubble.branches, key=len)[len(bubble.branches) // 2])
        if len(median_branch) == 0:
            continue

        #Bubble is TOO BIIG, will not correct it (maybe at the next iteration)
        if len(median_branch) > MAX_BUBBLE * 1.5:
            new_branches = [median_branch]
            long_branches += 1

        else:
            for branch in bubble.branches:
                incons_rate = abs(len(branch) - len(median_branch)) / len(median_branch)
                if incons_rate < 0.5:
                    if len(branch) == 0:
                        branch = "A"
                        #logger.debug("Zero branch")
                    new_branches.append(branch)

        if (abs(len(median_branch) - len(bubble.consensus)) > len(median_branch) // 2):
            bubble.consensus = median_branch

        if len(new_branches) > MAX_BRANCHES:
            new_branches = new_branches[:MAX_BRANCHES]

        new_bubbles.append(Bubble(bubble.contig_id, bubble.position))
        new_bubbles[-1].consensus = bubble.consensus
        new_bubbles[-1].branches = new_branches

    return new_bubbles, empty_bubbles, long_branches


def _is_solid_kmer(profile, position, err_mode):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = cfg.vals["err_modes"][err_mode]["solid_missmatch"]
    INS_RATE = cfg.vals["err_modes"][err_mode]["solid_indel"]
    SOLID_LEN = cfg.vals["solid_kmer_length"]

    for i in range(position, position + SOLID_LEN):
        if profile[i].coverage == 0:
            return False
        local_missmatch = (profile[i].num_missmatch +
                           profile[i].num_deletions) / profile[i].coverage
        local_ins = profile[i].num_inserts / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_simple_kmer(profile, position):
    """
    Checks if the kmer with center at the given position is simple
    """
    SIMPLE_LEN = cfg.vals["simple_kmer_length"]

    extended_len = SIMPLE_LEN * 2
    nucl_str = [p.nucl for p in profile[position - extended_len // 2 :
                                        position + extended_len // 2]]

    #single nucleotide homopolymers
    for i in range(extended_len // 2 - SIMPLE_LEN // 2,
                   extended_len // 2 + SIMPLE_LEN // 2 - 1):
        if nucl_str[i] == nucl_str[i + 1]:
            return False

    #dinucleotide homopolymers
    for shift in [0, 1]:
        for i in range(SIMPLE_LEN - shift - 1):
            pos = extended_len // 2 - SIMPLE_LEN + shift + i * 2
            if (nucl_str[pos : pos + 2] == nucl_str[pos + 2 : pos + 4]):
                return False

    #trinucleotide homopolymers
    #for shift in [0, 1, 2]:
    #    for i in xrange(SIMPLE_LEN - shift - 1):
    #        pos = shift + i * 3
    #        if (nucl_str[pos : pos + 3] == nucl_str[pos + 3 : pos + 6]):
    #            #logger.debug("tri" + "".join(nucl_str))
    #            return False

    return True


def _compute_profile(alignment, platform, genome_len):
    """
    Computes alignment profile
    """
    max_aln_err = cfg.vals["err_modes"][platform]["max_aln_error"]
    min_aln_len = cfg.vals["min_polish_aln_len"]
    aln_errors = []
    #filtered = 0
    profile = [ProfileInfo() for _ in range(genome_len)]
    for aln in alignment:
        if aln.err_rate > max_aln_err or len(aln.qry_seq) < min_aln_len:
            #filtered += 1
            continue

        aln_errors.append(aln.err_rate)

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in zip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-":
                prof_elem.num_inserts += 1
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.coverage += 1

                if qry_nuc == "-":
                    prof_elem.num_deletions += 1
                elif trg_nuc != qry_nuc:
                    prof_elem.num_missmatch += 1

            trg_pos += 1

    #logger.debug("Filtered: {0} out of {1}".format(filtered, len(alignment)))
    return profile, aln_errors


def _get_partition(profile, err_mode):
    """
    Partitions genome into sub-alignments at solid regions / simple kmers
    """
    #logger.debug("Partitioning genome")
    SOLID_LEN = cfg.vals["solid_kmer_length"]
    SIMPLE_LEN = cfg.vals["simple_kmer_length"]
    MAX_BUBBLE = cfg.vals["max_bubble_length"]

    solid_flags = [False for _ in range(len(profile))]
    prof_pos = 0
    while prof_pos < len(profile) - SOLID_LEN:
        if _is_solid_kmer(profile, prof_pos, err_mode):
            for i in range(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    partition = []
    prev_partition = SOLID_LEN

    long_bubbles = 0
    prof_pos = SOLID_LEN
    while prof_pos < len(profile) - SOLID_LEN:
        cur_partition = prof_pos + SIMPLE_LEN // 2
        landmark = (all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]) and
                    _is_simple_kmer(profile, cur_partition))

        if prof_pos - prev_partition > MAX_BUBBLE:
            long_bubbles += 1

        if landmark or prof_pos - prev_partition > MAX_BUBBLE:
            partition.append(cur_partition)
            prev_partition = cur_partition
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    #logger.debug("Partitioned into {0} segments".format(len(partition) + 1))
    #logger.debug("Long bubbles: {0}".format(long_bubbles))

    return partition, long_bubbles


def _get_bubble_seqs(alignment, platform, profile, partition, contig_info):
    """
    Given genome landmarks, forms bubble sequences
    """
    if not partition:
        return []

    #max_aln_err = cfg.vals["err_modes"][platform]["max_aln_error"]
    bubbles = []
    ext_partition = [0] + partition + [contig_info.length]
    for p_left, p_right in zip(ext_partition[:-1], ext_partition[1:]):
        bubbles.append(Bubble(contig_info.id, p_left))
        consensus = [p.nucl for p in profile[p_left : p_right]]
        bubbles[-1].consensus = "".join(consensus)

    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue

        bubble_id = bisect(partition, aln.trg_start % contig_info.length)
        next_bubble_start = ext_partition[bubble_id + 1]
        chromosome_start = (bubble_id == 0 and
                            not contig_info.type == "circular")
        chromosome_end = (aln.trg_end > partition[-1] and not
                          contig_info.type == "circular")

        branch_start = None
        first_segment = True
        trg_pos = aln.trg_start
        for i, trg_nuc in enumerate(aln.trg_seq):
            if trg_nuc == "-":
                continue
            if trg_pos >= contig_info.length:
                trg_pos -= contig_info.length

            if trg_pos >= next_bubble_start or trg_pos == 0:
                if not first_segment or chromosome_start:
                    branch_seq = fp.to_acgt(aln.qry_seq[branch_start : i].replace("-", ""))
                    bubbles[bubble_id].branches.append(branch_seq)

                first_segment = False
                bubble_id = bisect(partition, trg_pos)
                next_bubble_start = ext_partition[bubble_id + 1]
                branch_start = i

            trg_pos += 1

        if chromosome_end:
            branch_seq = fp.to_acgt(aln.qry_seq[branch_start:].replace("-", ""))
            bubbles[-1].branches.append(branch_seq)

    return bubbles
