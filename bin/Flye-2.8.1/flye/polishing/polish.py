#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs polishing binary in parallel and concatentes output
"""

from __future__ import absolute_import
from __future__ import division
import logging
import subprocess
import os
from collections import defaultdict

from flye.polishing.alignment import (make_alignment, get_contigs_info,
                                      merge_chunks, split_into_chunks)
from flye.utils.sam_parser import SynchronizedSamReader
from flye.polishing.bubbles import make_bubbles
import flye.utils.fasta_parser as fp
from flye.utils.utils import which
import flye.config.py_cfg as cfg
from flye.six import iteritems
from flye.six.moves import range


POLISH_BIN = "flye-modules"

logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polishing binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([POLISH_BIN, "polisher", "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def polish(contig_seqs, read_seqs, work_dir, num_iters, num_threads, error_mode,
           output_progress):
    """
    High-level polisher interface
    """
    logger_state = logger.disabled
    if not output_progress:
        logger.disabled = True

    subs_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][error_mode]["subs_matrix"])
    hopo_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][error_mode]["hopo_matrix"])
    stats_file = os.path.join(work_dir, "contigs_stats.txt")

    prev_assembly = contig_seqs
    contig_lengths = None
    coverage_stats = None
    for i in range(num_iters):
        logger.info("Polishing genome (%d/%d)", i + 1, num_iters)

        #split into 1Mb chunks to reduce RAM usage
        #slightly vary chunk size between iterations
        CHUNK_SIZE = 1000000 - (i % 2) * 100000
        chunks_file = os.path.join(work_dir, "chunks_{0}.fasta".format(i + 1))
        chunks = split_into_chunks(fp.read_sequence_dict(prev_assembly),
                                       CHUNK_SIZE)
        fp.write_fasta_dict(chunks, chunks_file)

        ####
        logger.info("Running minimap2")
        alignment_file = os.path.join(work_dir, "minimap_{0}.bam".format(i + 1))
        make_alignment(chunks_file, read_seqs, num_threads,
                       work_dir, error_mode, alignment_file,
                       reference_mode=True, sam_output=True)

        #####
        logger.info("Separating alignment into bubbles")
        contigs_info = get_contigs_info(chunks_file)
        bubbles_file = os.path.join(work_dir,
                                    "bubbles_{0}.fasta".format(i + 1))
        coverage_stats, mean_aln_error = \
            make_bubbles(alignment_file, contigs_info, chunks_file,
                         error_mode, num_threads,
                         bubbles_file)

        logger.info("Alignment error rate: %f", mean_aln_error)
        consensus_out = os.path.join(work_dir, "consensus_{0}.fasta".format(i + 1))
        polished_file = os.path.join(work_dir, "polished_{0}.fasta".format(i + 1))
        if os.path.getsize(bubbles_file) == 0:
            logger.info("No reads were aligned during polishing")
            if not output_progress:
                logger.disabled = logger_state
            open(stats_file, "w").write("#seq_name\tlength\tcoverage\n")
            open(polished_file, "w")
            return polished_file, stats_file

        #####
        logger.info("Correcting bubbles")
        _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                        consensus_out, num_threads, output_progress)
        polished_fasta, polished_lengths = _compose_sequence(consensus_out)
        merged_chunks = merge_chunks(polished_fasta)
        fp.write_fasta_dict(merged_chunks, polished_file)

        #Cleanup
        os.remove(chunks_file)
        os.remove(bubbles_file)
        os.remove(consensus_out)
        os.remove(alignment_file)

        contig_lengths = polished_lengths
        prev_assembly = polished_file

    #merge information from chunks
    contig_lengths = merge_chunks(contig_lengths, fold_function=sum)
    coverage_stats = merge_chunks(coverage_stats,
                                  fold_function=lambda l: sum(l) // len(l))

    with open(stats_file, "w") as f:
        f.write("#seq_name\tlength\tcoverage\n")
        for ctg_id in contig_lengths:
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    contig_lengths[ctg_id], coverage_stats[ctg_id]))

    if not output_progress:
        logger.disabled = logger_state

    return prev_assembly, stats_file


def generate_polished_edges(edges_file, gfa_file, polished_contigs, work_dir,
                            error_mode, num_threads):
    """
    Generate polished graph edges sequences by extracting them from
    polished contigs
    """
    logger.debug("Generating polished GFA")

    alignment_file = os.path.join(work_dir, "edges_aln.bam")
    polished_dict = fp.read_sequence_dict(polished_contigs)
    make_alignment(polished_contigs, [edges_file], num_threads,
                   work_dir, error_mode, alignment_file,
                   reference_mode=True, sam_output=True)
    aln_reader = SynchronizedSamReader(alignment_file,
                                       polished_dict,
                                       cfg.vals["max_read_coverage"])
    aln_by_edge = defaultdict(list)

    #getting one best alignment for each contig
    while not aln_reader.is_eof():
        _, ctg_aln = aln_reader.get_chunk()
        for aln in ctg_aln:
            aln_by_edge[aln.qry_id].append(aln)
    aln_reader.close()

    MIN_CONTAINMENT = 0.9
    updated_seqs = 0
    edges_dict = fp.read_sequence_dict(edges_file)
    for edge in edges_dict:
        if edge in aln_by_edge:
            main_aln = aln_by_edge[edge][0]
            map_start = main_aln.trg_start
            map_end = main_aln.trg_end
            for aln in aln_by_edge[edge]:
                if aln.trg_id == main_aln.trg_id and aln.trg_sign == main_aln.trg_sign:
                    map_start = min(map_start, aln.trg_start)
                    map_end = max(map_end, aln.trg_end)

            new_seq = polished_dict[main_aln.trg_id][map_start : map_end]
            if main_aln.qry_sign == "-":
                new_seq = fp.reverse_complement(new_seq)

            #print edge, main_aln.qry_len, len(new_seq), main_aln.qry_start, main_aln.qry_end
            if len(new_seq) / aln.qry_len > MIN_CONTAINMENT:
                edges_dict[edge] = new_seq
                updated_seqs += 1

    #writes fasta file with polished egdes
    #edges_polished = os.path.join(work_dir, "polished_edges.fasta")
    #fp.write_fasta_dict(edges_dict, edges_polished)

    #writes gfa file with polished edges
    with open(os.path.join(work_dir, "polished_edges.gfa"), "w") as gfa_polished, \
         open(gfa_file, "r") as gfa_in:
        for line in gfa_in:
            if line.startswith("S"):
                seq_id = line.split()[1]
                coverage_tag = line.split()[3]
                gfa_polished.write("S\t{0}\t{1}\t{2}\n"
                                    .format(seq_id, edges_dict[seq_id], coverage_tag))
            else:
                gfa_polished.write(line)

    logger.debug("%d sequences remained unpolished",
                 len(edges_dict) - updated_seqs)
    os.remove(alignment_file)


def filter_by_coverage(args, stats_in, contigs_in, stats_out, contigs_out):
    """
    Filters out contigs with low coverage
    """
    SUBASM_MIN_COVERAGE = 1
    HARD_MIN_COVERAGE = cfg.vals["hard_minimum_coverage"]
    RELATIVE_MIN_COVERAGE = cfg.vals["relative_minimum_coverage"]

    ctg_stats = {}
    sum_cov = 0
    sum_length = 0

    with open(stats_in, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            tokens = line.split("\t")
            ctg_id, ctg_len, ctg_cov = tokens[0], int(tokens[1]), int(tokens[2])
            ctg_stats[ctg_id] = (ctg_len, ctg_cov)
            sum_cov += ctg_cov * ctg_len
            sum_length += ctg_len

    mean_coverage = int(sum_cov / sum_length)
    coverage_threshold = None
    if args.read_type == "subasm":
        coverage_threshold = SUBASM_MIN_COVERAGE
    elif args.meta:
        coverage_threshold = HARD_MIN_COVERAGE
    else:
        coverage_threshold = int(round(mean_coverage /
                                       RELATIVE_MIN_COVERAGE))
        coverage_threshold = max(HARD_MIN_COVERAGE, coverage_threshold)
    logger.debug("Mean contig coverage: %d, selected threshold: %d",
                 mean_coverage, coverage_threshold)

    filtered_num = 0
    filtered_seq = 0
    good_fasta = {}
    for hdr, seq in fp.stream_sequence(contigs_in):
        if ctg_stats[hdr][1] >= coverage_threshold:
            good_fasta[hdr] = seq
        else:
            filtered_num += 1
            filtered_seq += ctg_stats[hdr][0]
    logger.debug("Filtered %d contigs of total length %d",
                 filtered_num, filtered_seq)

    fp.write_fasta_dict(good_fasta, contigs_out)
    with open(stats_out, "w") as f:
        f.write("#seq_name\tlength\tcoverage\n")
        for ctg_id in good_fasta:
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    ctg_stats[ctg_id][0], ctg_stats[ctg_id][1]))


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads, output_progress):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, "polisher", "--bubbles", bubbles_in, "--subs-mat", subs_matrix,
               "--hopo-mat", hopo_matrix, "--out", consensus_out,
               "--threads", str(num_threads)]
    if not output_progress:
        cmdline.append("--quiet")

    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def _compose_sequence(consensus_file):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    coverage = defaultdict(list)
    with open(consensus_file, "r") as f:
        header = True
        for line in f:
            if header:
                tokens = line.strip().split(" ")
                ctg_id = tokens[0][1:]
                ctg_pos = int(tokens[1])
                coverage[ctg_id].append(int(tokens[2]))
            else:
                consensuses[ctg_id].append((ctg_pos, line.strip()))
            header = not header

    polished_fasta = {}
    polished_stats = {}
    for ctg_id, seqs in iteritems(consensuses):
        sorted_seqs = [p[1] for p in sorted(seqs, key=lambda p: p[0])]
        concat_seq = "".join(sorted_seqs)
        #mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        polished_fasta[ctg_id] = concat_seq
        polished_stats[ctg_id] = len(concat_seq)

    return polished_fasta, polished_stats
