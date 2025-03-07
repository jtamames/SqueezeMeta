#!/usr/bin/env python
#Fernando Puente-Sánchez & Matthias Hoetzinger, 2023

import math
from argparse import ArgumentParser
from subprocess import run
from os import mkdir
from os.path import isdir
from sys import exit
from shutil import rmtree
from pandas import read_csv
from multiprocessing import Pool
from functools import partial


def main(args):
    inputVCF   = args.inputVCF
    outDir     = args.output
    minCount   = args.minCount
    nReps      = args.replicates
    processors = args.threads
    GFF        = args.gff
    inputFasta = args.fasta
    prefix     = args.prefix

    if(bool(GFF) ^ bool(inputFasta)): # XOR operator
       raise Exception('--gff and --inputFasta must either be provided together or not at all')


    # Derived params
    tempDir = f'{outDir}/temp'
    includedGFF = f'{tempDir}/selected.gff'

    nSamples = 0
    with open(inputVCF) as infile:
        for line in infile:
            if line.startswith('#CHROM'):
                chrom, pos, id_, ref, alt, qual, filter_, info, format_, *samples = line.strip().split('\t')
                nSamples = len(samples)
                break
    assert nSamples


    # Create output dir
    if not isdir(outDir):
        mkdir(outDir)
    if not isdir(tempDir):
        mkdir(tempDir)

    if GFF:
        # Prepare gff file
        print('Filtering gff file...')
        includedContigs = {line.replace('##contig=<ID=', '').split(',length=')[0] for line in open(inputVCF) if line.startswith('##contig=<ID=')}
        with open(GFF) as infile, open(includedGFF, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    continue
                contig = line.split('\t')[0]
                feature = line.split('\t')[2]
                if contig in includedContigs and feature == 'CDS':
                    outfile.write(line)


    # Remove low confidence sites, keep only high confidence poly and monomorphic loci
    print('Removing low confidence loci...')
    removed_mono, removed_poly, kept_mono, kept_poly = filter_vcf(0.99, 0.99, inputVCF, f'{tempDir}/filtered.vcf', verbose = True)

    if not kept_poly:
        print('Your VCF file contains no high confidence polymorphic sites, so we can not proceed. Stopping.')
        exit(0)


    # Run pogenom.pl to get the right number of good loci (high confidence + passing the min_found filter)
    print('Calculating number of good loci...')
    run(['pogenom.pl',
         '-vcf_file', f'{tempDir}/filtered.vcf',
         '--out', f'{tempDir}/POGENOM_prerun',
         '--min_count', str(minCount),
         '--subsample', str(minCount),
         '--min_found', str(nSamples),
         '--pi_only',
         '--genome_size', str(50000) # mock number here, since we don't care about the actual result
        ], capture_output = True)

    numLoci = 0
    with open(f'{tempDir}/POGENOM_prerun.intradiv.txt') as infile:
        for line in infile:
            if line.startswith('Sample'):
                continue
            numLoci = int(line.strip().split('\t')[4])
            break
    assert numLoci


    # Now remove the monomorphic loci for extra speed
    print('Removing monomorphic loci')
    removed_mono, removed_poly, kept_mono, kept_poly = filter_vcf(1, 0, f'{tempDir}/filtered.vcf', f'{tempDir}/filtered_nomono.vcf', verbose = True)

    # And do the real POGENOM run (have to do shell=True, otherwise it doesn't go through aminoacid stuff...
    msg = 'Running POGENOM for pi and gene-wide statistics' if GFF else 'Running POGENOM for pi calculation'
    command = ['pogenom.pl',
               '--vcf_file', f'{tempDir}/filtered_nomono.vcf',
               '--out', f'{outDir}/{prefix}',
               '--min_count', str(minCount),
               '--subsample', str(minCount),
               '--min_found', str(nSamples),
               '--genome_size', str(numLoci),
               '--fst_perm', str(nReps)
               ]
    if GFF:
        with open(f'{tempDir}/standard_genetic_code.txt', 'w') as outfile:
            outfile.write(GENCODE)
        command += ['--fasta_file', inputFasta, '--gff_file ', includedGFF, '--genetic_code_file', f'{tempDir}/standard_genetic_code.txt']
    
    print(msg)
    run(' '.join(command), shell = True, capture_output = True)

    # Now do replicates for FST
    if nSamples > 1:
        print(f'Running POGENOM in parallel for Fst calculation over {nReps} replicates')
        pool = Pool(processors)
        runrepP = partial(runRep, tempDir = tempDir, minCount = minCount, nSamples = nSamples, numLoci = numLoci)
        pool.map(runrepP, range(nReps))

        fstMinFoundAll = sum([read_csv(f'{tempDir}/{i}.min_foundA.fst.txt', sep='\t', index_col = 0) for i in range(nReps)]) / nReps
        fstMinFoundOne = sum([read_csv(f'{tempDir}/{i}.min_found1.fst.txt', sep='\t', index_col = 0) for i in range(nReps)]) / nReps

        fstMinFoundAll.round(4).to_csv(f'{outDir}/{prefix}.fst.min_found_{nSamples}_replicates_{nReps}.txt',
                                       sep = '\t', na_rep = 'NA')
        fstMinFoundOne.round(4).to_csv(f'{outDir}/{prefix}.fst.min_found_1_replicates_{nReps}.txt',
                                       sep = '\t', na_rep = 'NA')
    #yes I think Anders used min_found 1 in the Pogenom paper.. but I would recommend to set min_found to the number of metagenomes bc then you have the same loci analyzed for each metagenome



def filter_vcf(mono_prob_cutoff, poly_prob_cutoff, file_vcf, filtered_file_vcf, verbose = True):
    # Matthias Hoetzinger, 2023

    # filters vcf file based on QUAL
    # uses separate cutoffs for monomorphic and polymorphic sites that are given by the user
    # uncomment part in line 28 and lines 48 and 59 if removed sites should be written to files

    ############################################### USAGE ###############################################
    # python3 filter_vcf.py <mono_prob_cutoff> <poly_prob_cutoff> <file.vcf> <filtered_file.vcf>

    #### EXAMPLE 1 (keep only sites with >99% probability of beeing mono- and polymorphic, respectively)
    # python3 filter_vcf.py 0.99 0.99 file.vcf filtered_file.vcf

    #### EXAMPLE 2 (remove all monomporphic sites and keep all polymorphic)
    # python3 filter_vcf.py 1 0 file.vcf filtered_file.vcf
    #####################################################################################################


    #### calculate QUAL cutoffs from user input
    # QUAL = -10 * logP
    # where P is the probability that site is not polymorphic (probability that site is polymorphic is 1-P)
    mono_cutoff = -10*math.log10(float(mono_prob_cutoff))	#e.g. user input of 0.99 => site will be considered if QUAL < 0.043648
    poly_cutoff = -10*math.log10(1-float(poly_prob_cutoff))	#e.g. user input of 0.99 => site will be considered if QUAL > 20

    # counters for removed (and kept) monomorphic and polymorphic sites
    removed_mono = 0
    removed_poly = 0
    kept_mono    = 0
    kept_poly    = 0

    with open(file_vcf, 'r') as vcf_in, open(filtered_file_vcf, 'w') as filtered_vcf: #, open('removed_mono.txt', 'w') as rem_mono, open('removed_poly.txt', 'w') as rem_poly:
            
        for line in vcf_in:

            #only work on non-empty lines
            if line.strip():

                #write the header lines
                if line.startswith('#'):
                    filtered_vcf.write(line)

                #### work on monomorphic sites
                elif line.split('\t')[4] == '.':

                    #write to file if QUAL is smaller than cutoff for monomorphic sites
                    if float(line.split('\t')[5]) < mono_cutoff:
                        filtered_vcf.write(line)
                        kept_mono += 1

                    #write to removed_mono.txt if QUAL isn't smaller and count
                    else:
                        #rem_mono.write(line)
                        removed_mono += 1

                #### work on polymorphic sites (all that are not monomorphic)
                else:
                    #write to file if QUAL is higher than cutoff for polymorphic sites
                    if float(line.split('\t')[5]) > poly_cutoff:
                        filtered_vcf.write(line)
                        kept_poly += 1
                    #write to removed_ploy.txt if QUAL isn't higher and count
                    else:
                        #rem_poly.write(line)
                        removed_poly += 1
    if verbose:
        print('#'*80+'\nRemoved '+str(removed_mono)+' monomorphic and '+str(removed_poly)+' polymorphic sites.')
        print('#'*80+'\nKept '+str(kept_mono)+' monomorphic and '+str(kept_poly)+' polymorphic sites.')
        print('Wrote filtered vcf file to \"'+filtered_file_vcf+'\".\n'+'#'*80)

    return removed_mono, removed_poly, kept_mono, kept_poly




def runRep(i, tempDir, minCount, nSamples, numLoci):
    run(['pogenom.pl',
         '--vcf_file', f'{tempDir}/filtered_nomono.vcf',
         '--out', f'{tempDir}/{i}.min_foundA',
         '--min_count', str(minCount),
         '--subsample', str(minCount),
         '--min_found', str(nSamples),
         '--genome_size', str(numLoci)
         ], capture_output = True)
    run(['pogenom.pl',
         '--vcf_file', f'{tempDir}/filtered_nomono.vcf',
         '--out', f'{tempDir}/{i}.min_found1',
         '--min_count', str(minCount),
         '--subsample', str(minCount),
         '--min_found', '1',
         '--genome_size', str(numLoci)
        ], capture_output = True)




def parse_args():
    parser = ArgumentParser(description='Run pogenom.pl after a successful run of the Input_POGENOM pipeline')
    parser.add_argument('inputVCF', type = str,
                        help = 'Path to the VCF file produced by Input_POGENOM.sh')
    parser.add_argument('-o', '--output', type = str, default = 'pogenom_out',
                        help = 'Output directory (will be created if inexistent)')
    parser.add_argument('-m', '--minCount', type = int, default = 10,
                        help = 'Minimum coverage (integer) for a locus to be included')
    parser.add_argument('-r', '--replicates', type = int, default = 100,
                        help = 'Replicates/permutations for Fst/randomised gene-wise Fst')
    parser.add_argument('-g', '--gff', type = str,
                        help = 'Path to a master GFF file including gene calls from the contigs in the input genome (required for gene-wise statistics)')
    parser.add_argument('-f', '--fasta', type = str,
                        help = 'Path to the fasta file used to run Input_POGENOM.sh (required for gene-wise statistics)')
    parser.add_argument('-p', '--prefix', type = str, default = 'POGENOM',
                       help = 'Output file prefix')
    parser.add_argument('-t', '--threads', type = int, default = 12,
                        help = 'Processors used for Fst replicate calculation')
    return parser.parse_args()


GENCODE = """TTT\tF\tPhe
TTC\tF\tPhe
TTA\tL\tLeu
TTG\tL\tLeu
CTT\tL\tLeu
CTC\tL\tLeu
CTA\tL\tLeu
CTG\tL\tLeu
ATT\tI\tIle
ATC\tI\tIle
ATA\tI\tIle
ATG\tM\tMet
GTT\tV\tVal
GTC\tV\tVal
GTA\tV\tVal
GTG\tV\tVal
TCT\tS\tSer
TCC\tS\tSer
TCA\tS\tSer
TCG\tS\tSer
CCT\tP\tPro
CCC\tP\tPro
CCA\tP\tPro
CCG\tP\tPro
ACT\tT\tThr
ACC\tT\tThr
ACA\tT\tThr
ACG\tT\tThr
GCT\tA\tAla
GCC\tA\tAla
GCA\tA\tAla
GCG\tA\tAla
TAT\tY\tTyr
TAC\tY\tTyr
TAA\t*\tTer
TAG\t*\tTer
CAT\tH\tHis
CAC\tH\tHis
CAA\tQ\tGln
CAG\tQ\tGln
AAT\tN\tAsn
AAC\tN\tAsn
AAA\tK\tLys
AAG\tK\tLys
GAT\tD\tAsp
GAC\tD\tAsp
GAA\tE\tGlu
GAG\tE\tGlu
TGT\tC\tCys
TGC\tC\tCys
TGA\t*\tTer
TGG\tW\tTrp
CGT\tR\tArg
CGC\tR\tArg
CGA\tR\tArg
CGG\tR\tArg
AGT\tS\tSer
AGC\tS\tSer
AGA\tR\tArg
AGG\tR\tArg
GGT\tG\tGly
GGC\tG\tGly
GGA\tG\tGly
GGG\tG\tGly""".replace('\n', '\r\n') # this produces the same md5sum as the original file from POGENOM
                                     #  which had DOS line terminators

if __name__ == '__main__':
    main(parse_args())
