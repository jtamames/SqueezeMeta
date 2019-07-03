#!/usr/bin/python3

"""
Part of SqueezeMeta distribution. 25/03/2018 Original version, (c) Fernando Puente-Sánchez, CNB-CSIC

Create the files required for loading a SqueezeMeta project into the web interface (https://github.com/jtamames/SqueezeMdb).

USAGE: make-SqueezeMdb-files.py <PROJECT_NAME> <OUTPUT_DIRECTORY> 
"""


from os.path import abspath, dirname, realpath
from os import mkdir, system
from sys import exit, path
import argparse

utils_home = abspath(dirname(realpath(__file__)))
path.append('{}/../lib/'.format(utils_home))
from utils import parse_conf_file, write_orf_seqs

def main(args):
    ### Get result files paths from SqueezeMeta_conf.pl
    perlVars = parse_conf_file(args.project_path)

    ### Create output dir.
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise
        else:
            print('\nThe directory {} already exists. Please remove it or use a different output name.\n'.format(args.output_dir))
            exit(1)

    
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
    def new2old(str_, orftable=False):
        """Replace 1.0 headers with old headers, so we don't have to modify and re-deploy SQMdb"""
        if orftable:
            str_ = str_.replace('Coverage', 'COVERAGE').replace('Raw read count', 'RAW READ COUNT').replace('Raw base count', 'RAW BASE COUNT')
        else:
            str_ = str_.replace('Strain het', 'Strain Het').replace('Raw read count', 'Raw')
        return str_

    allORFs = []
    newFields =  ['ORF ID', 'Contig ID', 'Length AA', 'GC perc', 'Gene name', 'Tax', 'KEGG ID', 'KEGGFUN', 'KEGGPATH', 'COG ID', 'COGFUN', 'COGPATH', 'PFAM']
    goodFields = ['ORF', 'CONTIG ID', 'LENGTH AA', 'GC perc', 'GENNAME', 'TAX ORF', 'KEGG ID', 'KEGGFUN', 'KEGGPATH', 'COG ID', 'COGFUN', 'COGPATH', 'PFAM']
    with open(perlVars['$mergedfile']) as infile, open('{}/genes.tsv'.format(args.output_dir), 'w') as outfile:
        outfile.write(infile.readline())
        header = infile.readline().strip().split('\t')
        newFields.extend([f for f in header if f.startswith('TPM ') or f.startswith('Coverage ') or f.startswith('Raw read ') or f.startswith('Raw base ')])
        goodFields.extend([new2old(f, True) for f in header if f.startswith('TPM ') or f.startswith('Coverage ') or f.startswith('Raw read ') or f.startswith('Raw base ')])
        outfile.write('\t'.join(goodFields) + '\n')
        idx =  {f: i for i,f in enumerate(header) if f in newFields}
        for line in infile:
            line = line.strip().split('\t')
            if line[2] == 'CDS':
                allORFs.append(line[0])
                outfile.write('{}\n'.format('\t'.join([line[idx[f]] for f in newFields])))


    ### Create contigtable.
    with open(perlVars['$contigtable']) as infile, open('{}/contigs.tsv'.format(args.output_dir), 'w') as outfile:
        outfile.write(infile.readline())
        outfile.write(new2old(infile.readline())) # adapt header
        [outfile.write(line) for line in infile]
    

    ### Create bintable.
    if not int(perlVars['$nobins']):
        with open(perlVars['$bintable']) as infile, open('{}/bins.tsv'.format(args.output_dir), 'w') as outfile:
            outfile.write(infile.readline())
            outfile.write(new2old(infile.readline())) # adapt header
            [outfile.write(line) for line in infile]


    ### Create sequences file.
    aafile = perlVars['$aafile']
    fna_blastx = perlVars['$fna_blastx'] if int(perlVars['$doublepass']) else None
    outname = '{}/sequences.tsv'.format(args.output_dir)
    write_orf_seqs(allORFs, aafile, fna_blastx, None, outname)


def parse_args():
    parser = argparse.ArgumentParser(description='Create the files required for loading a SqueezeMeta project into the web interface', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')

    return parser.parse_args()


if __name__ == '__main__':
    main(parse_args())


