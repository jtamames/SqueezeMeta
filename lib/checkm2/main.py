#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import sys
import logging
import shutil
import tempfile
import pandas as pd
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import gzip
import tarfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from checkm2 import version

__author__ = "Alex Chklovski"
__version__ = version.__version__
__maintainer__ = "Alex Chklovski"
__email__ = "chklovski near gmail.com"
__status__ = "Development"

from checkm2.defaultValues import DefaultValues
from checkm2.versionControl import VersionControl
from checkm2 import fileManager
from checkm2 import predictQuality



def main():

    num_threads = 1


    class ChangeTempAction(argparse.Action):
        def __call__(self, parser, namespace, newtmpdir, option_string=None):
            if os.path.isdir(newtmpdir):
                tempfile.tempdir = newtmpdir
            else:
                raise argparse.ArgumentTypeError(
                    'The value of %s must be a valid directory' % option_string)


    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--lowmem', help='Low memory mode. Reduces DIAMOND blocksize to significantly reduce RAM usage at the expense of longer runtime', action="store_true", default=False)


    parser = argparse.ArgumentParser(parents=[parent_parser])
    subparsers = parser.add_subparsers(title="Sub-commands", dest='subparser_name', parser_class=argparse.ArgumentParser)
    subparser_name_to_parser = {}

    def new_subparser(subparsers, parser_name, parser_description):
        subpar = subparsers.add_parser(parser_name,
                                       description=parser_description,
                                       help=parser_description,
                                       formatter_class=RawTextHelpFormatter,
                                       parents=[parent_parser])
        subparser_name_to_parser[parser_name] = subpar
        return subpar


    predict_description = 'Predict the completeness and contamination of genome bins in a folder. Example usage: \n\n' \
                           '\tcheckm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder>'
    testrun_description = 'Runs Checkm2 on internal test genomes to ensure it runs without errors. Example usage: \n\n' \
                            '\t checkm2 testrun --threads 10'
    download_description = 'Download/set up required diamond database for CheckM2. Example usage: \n\n ' \
                           '\tcheckm2 database --download (downloads database into /home/user/databases)\n ' \
                           '\tcheckm2 database --download --path /path/to/custom_location  (downloads database into specified folder)\n ' \
                           '\tcheckm2 database --setdblocation /path/to/downloaded_database_file (uses specified database file as DB) \n\n ' \
                           'Alternatively, add an existing DIAMOND DB file to path: "export CHECKM2DB=/path/to/database/database.dmnd"\n\n'

    predict_parser = new_subparser(subparsers, 'predict', predict_description)

    predict_arguments = predict_parser.add_argument_group('required arguments')

    predict_arguments.add_argument('--input', '-i', help="Path to folder containing MAGs or list of MAGS to be analyzed", required=True, nargs='+')
    predict_arguments.add_argument('--output-directory', '--output_directory', '-o', help="Path output to folder", required=True)


    predict_arguments = predict_parser.add_argument_group('additional arguments')

    predict_arguments.add_argument('--general', action='store_true',
                                   help='Force the use of the general quality prediction model (gradient boost)', default=False)
    predict_arguments.add_argument('--specific', action='store_true',
                                   help='Force the use of the specific quality prediction model (neural network)', default=False)
    predict_arguments.add_argument('--allmodels', action='store_true',
                                   help='Output quality prediction for both models for each genome.', default=False)
    predict_arguments.add_argument('--genes', action='store_true',
                                   help='Treat input files as protein files. [Default: False]', default=False)
    predict_arguments.add_argument('-x', '--extension',
                                   help='Extension of input files. [Default: .fna]', default='.fna')
    predict_arguments.add_argument('--tmpdir', action=ChangeTempAction,
                                   help="specify an alternative directory for temporary files")

    predict_arguments.add_argument('--force', action='store_true', help='overwrite output directory [default: do not overwrite]', default=False)
    predict_arguments.add_argument('--resume', action='store_true', help='Reuse Prodigal and DIAMOND results found in output directory [default: not set]', default=False)
    predict_arguments.add_argument('--threads', '-t', type=int, metavar='num_threads', help='number of CPUS to use [default: %i]' % num_threads, default=num_threads)
    predict_arguments.add_argument('--stdout', action='store_true', help='Print results to stdout [default: write to file]', default=False)
    predict_arguments.add_argument('--remove_intermediates', action='store_true', help="Remove all intermediate files (protein files, diamond output) [default: don't]", default=False)
    predict_arguments.add_argument('--ttable', type=int, metavar='ttable', help="Provide a specific progidal translation table for bins [default: automatically determine either 11 or 4]", default=None)
    predict_arguments.add_argument('--database_path', help="Provide a location for the CheckM2 database for a given predict run [default: use either internal path set via <checkm2 database> or CHECKM2DB environmental variable]", default=None)

    predict_arguments.add_argument('--dbg_cos', action='store_true', help="DEBUG: write cosine similarity values to file [default: don't]", default=False)
    predict_arguments.add_argument('--dbg_vectors', action='store_true', help="DEBUG: dump pickled feature vectors to file [default: don't]", default=False)


    test_parser = new_subparser(subparsers, 'testrun', testrun_description)
    test_parser.add_argument('--threads', '-t', type=int, metavar='num_threads', help='number of CPUS to use [default: %i]' % num_threads, default=num_threads)
    test_parser.add_argument('--database_path', help="Provide a location for the CheckM2 database for a given predict run [default: use either internal path set via <checkm2 database> or CHECKM2DB environmental variable]", default=None)

    download_parser = new_subparser(subparsers, 'database', download_description)
    action  = download_parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--download', help="Download DIAMOND database. By default installs into [{}]".format(DefaultValues.DEFAULT_DB_INSTALL_LOCATION), action='store_true')
    action.add_argument('--setdblocation', help="Point CheckM2 to the DIAMOND database location if already downloaded.")
    action.add_argument('--current',  action='store_true', help="Print where current database is installed.")
    download_parser.add_argument('--path', help='Custom path for downloading and installing database file.', default=DefaultValues.DEFAULT_DB_INSTALL_LOCATION)
    download_parser.add_argument('--no_write_json_db',  help='Do NOT attempt to write database path to internal JSON file [useful if install directory is not writable]', action='store_true', default=False)



    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == 'help'):
        print('                ...::: CheckM2 v' + __version__ + ' :::...''')
        print('\n  General usage:')
        print('    predict         -> %s' % predict_description)
        print('    testrun         -> %s' % testrun_description)
        print('    database        -> %s' % 'Download and set up required CheckM2 DIAMOND database for annotation')

        print('\n  Use checkm2 <command> -h for command-specific help.\n')
        sys.exit(0)

    else:
        args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO


    logging.basicConfig(level=loglevel, format='[%(asctime)s] %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    def validate_args_model_choice(args):
        if (args.specific and args.general) or (args.specific and args.allmodels) or (args.general and args.allmodels):
            logging.error("Only one of --general --specific --allmodels can be specified.")
            sys.exit(1)
        if args.specific:
            return 'specific'
        elif args.general:
            return 'general'
        if args.allmodels:
            return 'both'
        else:
            return 'auto'

    if args.subparser_name == 'predict':

        #check if folder is empty and force remove it if necessary
        if not args.resume:
            fileManager.check_empty_dir(args.output_directory, args.force)
        else:
            fileManager.check_if_dir_exists(args.output_directory)

        #check custom translation tables

        allowable_ttables = [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33]

        if args.ttable is not None:
            if int(args.ttable) not in allowable_ttables:
                logging.error('Translation table {} is not valid'.format(args.ttable))
                sys.exit(1)


        logging.root.handlers = []
        logging.basicConfig(level=loglevel, format='[%(asctime)s] %(levelname)s: %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(os.path.join(os.path.abspath(args.output_directory), "checkm2.log")),
                                      logging.StreamHandler()])

        mode = validate_args_model_choice(args)


        logging.info('Running CheckM2 version {}'.format(__version__))

        # check if db location was provided
        tempDBpath = None
        if args.database_path is not None:
            logging.info(f'Custom database path provided for predict run. Checking database at {args.database_path}...')
            if VersionControl().checksum_version_validate_DIAMOND(args.database_path):
                tempDBpath = args.database_path

        logging.info("Running quality prediction workflow with {} threads.".format(args.threads))
        if len(args.input) == 1 and os.path.isdir(args.input[0]):
                predictor = predictQuality.Predictor(args.input[0], args.output_directory, args.extension, args.threads,
                                                     args.lowmem, tempDBpath)
                
                predictor.prediction_wf(args.genes, mode, args.dbg_cos, args.dbg_vectors, args.stdout,
                                        args.resume, args.remove_intermediates, args.ttable)
        else:
            if args.genes:
                bin_extension = 'faa'
            else:
                bin_extension = 'fna'
            # make folder to copy disparate bins into, then point to it as input folder
            # make them all one extension to make life easier
            # fileManager.make_sure_path_exists(os.path.join(args.output_directory, 'input_bins'))
            bin_temporary_dir = tempfile.TemporaryDirectory()
            for bin in args.input:
                if os.stat(bin).st_size == 0:
                    logging.warning("Skipping bin {} as it has a size of 0 bytes.".format(bin))
                elif os.path.isdir(bin):
                    continue
                elif bin.endswith('.gz'):
                    if tarfile.is_tarfile(bin):
                        logging.warning('Skipping file {} as tar archives are not supported'.format(bin))
                    else:
                        basename = os.path.splitext(os.path.basename(bin))[0]
                        #remove the .x.gz extension to have uniform labelling
                        if len(bin.split('.')) > 2:
                            pre_gz_prefix = '.{}'.format(bin.split('.')[-2])
                            basename = basename.split(pre_gz_prefix)[0]
                        with gzip.open(bin, 'rb') as f_in:
                            with open(os.path.join(bin_temporary_dir.name, '{}.{}'.format(basename, bin_extension)), 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                else:
                    shutil.copyfile(bin, os.path.join(bin_temporary_dir.name, '{}.{}'.format(os.path.splitext(os.path.basename(bin))[0], bin_extension)))
            predictor = predictQuality.Predictor(bin_temporary_dir.name, args.output_directory, bin_extension, args.threads,
                                                 args.lowmem, tempDBpath)
            predictor.prediction_wf(args.genes, mode, args.dbg_cos, args.dbg_vectors,
                                    args.stdout, args.resume, args.remove_intermediates, args.ttable)
            bin_temporary_dir.cleanup()

    elif args.subparser_name == 'testrun':
        logging.info("Test run: Running quality prediction workflow on test genomes with {} threads.".format(args.threads))
        logging.info('Running checksum on test genomes.')
        if VersionControl().checksum_test_genomes():
            logging.info('Checksum successful.')
        else:
            logging.error('Checksum unsuccessful. Please re-download test genomes or ensure CheckM2 version is correct.')
            sys.exit(1)

        # check if db location was provided
        tempDBpath = None
        if args.database_path is not None:
            logging.info(f'Custom database path provided for predict run. Checking database at {args.database_path}...')
            if VersionControl().checksum_version_validate_DIAMOND(args.database_path):
                tempDBpath = args.database_path

        with tempfile.TemporaryDirectory() as temp_out_dir:
            predictor = predictQuality.Predictor(DefaultValues.TESTRUN_GENOMES, temp_out_dir, '.tst', args.threads,
                                                 args.lowmem, tempDBpath)
            predictor.prediction_wf(False, 'auto', False, False, False)

            results = pd.read_csv(os.path.join(temp_out_dir, 'quality_report.tsv'), sep='\t')
            if len(results) > 1:
                logging.info('Test run successful! See README for details. Results:')
                print(results.iloc[:, :5].to_string(index=False, float_format=lambda x: '%.2f' % x))
            else:
                logging.error("Couldn't verify test run results.")


    elif args.subparser_name == 'database':

        if args.setdblocation:
            if not args.no_write_json_db:
                fileManager.DiamondDB().set_DB_location(args.setdblocation)
            else:
                logging.error('Cannot set db location when no_write_json_db flag is set.')
        elif args.download:
            fileManager.DiamondDB().download_database(args.path, args.no_write_json_db)
        elif args.current:
            loc = fileManager.DiamondDB().get_DB_location()
            logging.info(str(loc))

    else:
        raise Exception("Programming error")


if __name__ == '__main__':
    main()
