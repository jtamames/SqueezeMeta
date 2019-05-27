#!/usr/bin/python3

"""
Part of SqueezeMeta distribution. 25/03/2018 Original version, (c) Fernando Puente-Sánchez, CNB-CSIC

Create the files required for loading a SqueezeMeta project into the web interface (https://github.com/jtamames/SqueezeMdb).

USAGE: make-SqueezeMdb-files.py <PROJECT_NAME> <OUTPUT_DIRECTORY> 
"""

from os.path import abspath, dirname, realpath
from os import mkdir, system
from sys import path
import argparse

utils_home = abspath(dirname(realpath(__file__)))
path.append('{}/../lib/'.format(utils_home))
from utils import parse_conf_file

def main(args):
    ### Get result files paths from SqueezeMeta_conf.pl
    perlVars = parse_conf_file(args.project_path)

    ### Create output dir.
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise
    
    ### Create samples file.
    with open(perlVars['$mappingfile']) as infile, open('{}/samples.tsv'.format(args.output_dir), 'w') as outfile:
        outfile.write('Sample\ttest\n')
        addedSamples = set()
        for line in infile:
            sample = line.split('\t')[0].strip() # There shouldn't be trailing spaces though...
            if sample not in addedSamples:
                addedSamples.add(sample)
                outfile.write('{}\t{}\n'.format(sample, len(addedSamples)))


    ### Create orftable.
    allORFs = []
    goodFields = ['ORF', 'CONTIG ID', 'LENGTH AA', 'GC perc', 'GENNAME', 'TAX ORF', 'KEGG ID', 'KEGGFUN', 'KEGGPATH', 'COG ID', 'COGFUN', 'COGPATH', 'PFAM']
    with open(perlVars['$mergedfile']) as infile, open('{}/genes.tsv'.format(args.output_dir), 'w') as outfile:
        outfile.write(infile.readline())
        header = infile.readline().strip().split('\t')
        goodFields.extend([f for f in header if f.startswith('TPM ') or f.startswith('COVERAGE ') or f.startswith('RAW READ ') or f.startswith('RAW BASE ')])
        outfile.write('\t'.join(goodFields) + '\n')
        idx =  {f: i for i,f in enumerate(header) if f in goodFields}
        for line in infile:
            line = line.strip().split('\t')
            if line[2] == 'CDS':
                allORFs.append(line[0])
                outfile.write('{}\n'.format('\t'.join([line[idx[f]] for f in goodFields])))


    ### Create contigtable.
    system('cp {} {}/contigs.tsv'.format(perlVars['$contigtable'], args.output_dir))
    

    ### Create bintable.
    if not int(perlVars['$nobins']):
        system('cp {} {}/bins.tsv'.format(perlVars['$bintable'], args.output_dir))


    ### Create sequences file.
    # Load prodigal results.
    ORFseq = parse_fasta(perlVars['$aafile'])
    # Load blastx results if required.
    if int(perlVars['$doublepass']):
        ORFseq.update(parse_fasta(perlVars['$fna_blastx']))
    # Write results.
    with open('{}/sequences.tsv'.format(args.output_dir), 'w') as outfile:
        outfile.write('ORF\tAASEQ\n')
        for ORF in allORFs:
            outfile.write('{}\t{}\n'.format(ORF, ORFseq[ORF]))


def parse_fasta(fasta):
    """
    Parse a fasta file into a dictionary {header: sequence}
    """
    res = {}
    with open(fasta) as infile:
        header = ''
        seq = ''
        for line in infile:
            if line.startswith('>'):
                if header:
                    res[header] = seq
                    seq = ''
                header = line.strip().lstrip('>').split(' ')[0]
            else:
                seq += line.strip()
        res[header] = seq
    return res


def parse_args():
    parser = argparse.ArgumentParser(description='Create the files required for loading a SqueezeMeta project into the web interface', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')

    return parser.parse_args()


if __name__ == '__main__':
    main(parse_args())


