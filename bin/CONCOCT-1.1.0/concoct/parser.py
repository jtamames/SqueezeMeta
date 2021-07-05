from __future__ import print_function
import os
import sys
from random import randint
from argparse import ArgumentParser, ArgumentTypeError

def set_random_state(seed):
    ERROR="'{0}' should be converatable to integer".format(seed)
    try:
        seed = int(seed)
        if seed < 0:
            raise ArgumentTypeError("'" + seed + "' should be >= 0")
        elif seed == 0:
            seed = randint(2,10000)
        return seed
    except ValueError as e:
        raise ArgumentTypeError(ERROR)

def get_version():
  from concoct import __version__
  return '%(prog)s {version}'.format(version=__version__)

def arguments():
    parser = ArgumentParser()

    #Input files
    parser.add_argument('--coverage_file',
        help=("specify the coverage file, containing a table where each row "
              "correspond to a contig, and each column correspond to a sample. "
              "The values are the average coverage for this contig in that sample. "
              "All values are separated with tabs."))
    parser.add_argument('--composition_file',
        help=("specify the composition file, containing sequences in fasta format. "
              "It is named the composition file since it is used to calculate the "
              "kmer composition (the genomic signature) of each contig."))

    #Handle cluster number parsing
    parser.add_argument('-c', '--clusters', default=400, type=int,
      help='specify maximal number of clusters for VGMM, default 400.')
    #Kmer length, kmer count threshold and read length
    parser.add_argument('-k','--kmer_length', type=int, default=4,
        help='specify kmer length, default 4.')

    parser.add_argument('-t','--threads', type=int, default=1,
                    help='Number of threads to use')

    parser.add_argument('-l','--length_threshold', type=int, default=1000,
        help=("specify the sequence length threshold, contigs shorter than this "
              "value will not be included. Defaults to 1000."))
    parser.add_argument('-r','--read_length', type=int, default=100,
        help='specify read length for coverage, default 100')
    #Joined PCA
    parser.add_argument('--total_percentage_pca', default=90, type=int,
                        help=('The percentage of variance explained'
                              ' by the principal components for the'
                              ' combined data.'))
    #Output
    parser.add_argument('-b', '--basename', default=os.curdir,
      help=("Specify the basename for files or directory where output"
            "will be placed. Path to existing directory or basename"
            "with a trailing '/' will be interpreted as a directory."
            "If not provided, current directory will be used."))
    parser.add_argument('-s','--seed',type=set_random_state, default=set_random_state(1),
      help=('Specify an integer to use as seed for clustering. '
            '0 gives a random seed, 1 is the default seed and '
            'any other positive integer can be used. Other values '
            'give ArgumentTypeError.'))
    parser.add_argument('-i','--iterations',type=int, default=500,
      help=('Specify maximum number of iterations for the VBGMM. '
            'Default value is 500'))
    parser.add_argument('--no_cov_normalization', default=False, action="store_true",
      help=("By default the coverage is normalized with regards to samples, "
            "then normalized with regards of contigs and finally log transformed. "
            "By setting this flag you skip the normalization and only do log "
            "transorm of the coverage."))
    parser.add_argument('--no_total_coverage', default=False, action="store_true",
      help=("By default, the total coverage is added as a new column in the coverage "
            "data matrix, independently of coverage normalization but previous to "
            "log transformation. Use this tag to escape this behaviour."))
    parser.add_argument('--no_original_data', default=False, action="store_true",
      help=("By default the original data is saved to disk. For big datasets, "
            "especially when a large k is used for compositional data, this file can become "
            "very large. Use this tag if you don't want to save the original data."))
    parser.add_argument('-o','--converge_out', default=False, action="store_true",
      help=('Write convergence info to files.'))

    parser.add_argument('-d','--debug', default=False, action="store_true",
      help=('Debug parameters. '))
    parser.add_argument('-v','--version', action='version',
      version=get_version())

    args  = parser.parse_args()

    if args.debug:
        print(args, file=sys.stderr)
        sys.exit(0)
    # This can be changed to an or case when input of either case is supported individually
    if not (args.coverage_file or args.composition_file):
        parser.error("No input data supplied, add file(s) using --coverage_file <cov_file> and/or "
                     "--composition_file <comp_file>")

    return args
