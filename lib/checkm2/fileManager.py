import os
import errno
import sys
import logging
import shutil
import requests
import json
import gzip
import tempfile

from checkm2 import versionControl
from checkm2.defaultValues import DefaultValues
from checkm2 import zenodo_backpack

class DiamondDB:
    def __init__(self):
        if DefaultValues.DB_VAR in os.environ:
            self.DATABASE_DIR = os.environ[DefaultValues.DB_VAR]
            #Check if it's still there and if not, unset variable
            if not os.path.exists(self.DATABASE_DIR):
              logging.warning('Database not found using the environmental variable: {}. Please fix your $PATH. Using internal database path instead.'.format(DefaultValues.DB_VAR))

              diamond_definition = self.__get_db_file()

              if diamond_definition['DBPATH'] == 'Not Set':
                self.DATABASE_DIR = 'Not Set'
              else:
                self.DATABASE_DIR = diamond_definition['DBPATH']
        else:
            diamond_definition = self.__get_db_file()

            if diamond_definition['DBPATH'] == 'Not Set':
                self.DATABASE_DIR = 'Not Set'
            else:
                self.DATABASE_DIR = diamond_definition['DBPATH']


    def __get_db_file(self):
        diamond_location = DefaultValues.DB_LOCATION_DEFINITION
        try:
            with open(diamond_location) as f:
                diamond_definition = json.load(f)
                return diamond_definition
        except:
            logging.warning('Could not open DIAMOND location definition file. Please set up database or provide database path. Exiting.')
            sys.exit(1)


    def get_DB_location(self):
        if self.DATABASE_DIR == 'Not Set':
            logging.error('DIAMOND database not found. Please download database using <checkm2 database --download>')
            sys.exit(1)
        else:
            return self.DATABASE_DIR

    def set_DB_location(self, provided_location):
        logging.info('Checking provided DIAMOND database location')
        if versionControl.VersionControl().checksum_version_validate_DIAMOND(provided_location):
            #great - let's set it
            diamond_definition = self.__get_db_file()

            diamond_definition['DBPATH'] = os.path.abspath(provided_location)
            with open(DefaultValues.DB_LOCATION_DEFINITION, 'w') as dd:
                json.dump(diamond_definition, dd)

            logging.info("Database check successful! Database path successfully added.")

        else:
            logging.error("Checksum in CheckM2 reference doesn't match provided database. Please check your files.")
            sys.exit(1)


    def download_database(self, download_location, no_write_json_db):

        '''Uses a DOI link to automatically download, unpack and verify from zenodo.org'''

        logging.info("Command: Download database. Checking internal path information.")
        if no_write_json_db:
            logging.info('The --no_write_json_db flag was set. Will download database but NOT modify internal database path information. You will need to provide database path in the future.')

        diamond_location = DefaultValues.DB_LOCATION_DEFINITION
        
        try:
            with open(diamond_location) as f:
                diamond_definition = json.load(f)
        except:
            logging.warning('Could not open internal DIAMOND location definition file.')
            if not no_write_json_db: 
                logging.info('Writing new DIAMOND location definition file...')
                x = {"Type": "DIAMONDDB", "DBPATH": "Not Set"}
                with open(diamond_location, 'w') as dl:
                    json.dump(x, dl)
                try:
                    with open(diamond_location) as f:
                        diamond_definition = json.load(f)
                except Exception as e:
                    logging.error('Could not create new file: {}'.format(e))
                    sys.exit(1)
        if diamond_definition['DBPATH'] != 'Not Set':
            logging.warning('DIAMOND database found at {}.'.format(diamond_definition['DBPATH']))
            if not no_write_json_db:
                logging.info('Overwriting previous database')
            
        
        make_sure_path_exists(os.path.join(download_location, 'CheckM2_database'))
        
        backpack_downloader = zenodo_backpack.zenodo_backpack_downloader('INFO')
        highest_compatible_version, DOI = versionControl.VersionControl().return_highest_compatible_DB_version()

        
        diamond_loc_final = os.path.join(download_location, 'CheckM2_database', 'uniref100.KO.1.dmnd')

        if download_location is not None:
            #check we have writing permission
            try:
                os.makedirs(download_location, exist_ok=True)
                with tempfile.TemporaryDirectory(dir=download_location):
                    pass
            except OSError:
                logging.error("You do not appear to have permission to write to {}. Please choose a different directory"
                              .format(download_location))
                sys.exit(1)
            
            backpack_downloader.download_and_extract(download_location, DOI, progress_bar=True, no_check_version=False)
            
            diamond_definition['DBPATH'] = os.path.abspath(diamond_loc_final)
            if not no_write_json_db:
                with open(diamond_location, 'w') as dd:
                    json.dump(diamond_definition, dd)
            

        else:
            logging.info('Failed to determine download location')
            sys.exit(1)

        #do checksum

        if no_write_json_db:
            location = os.path.abspath(diamond_loc_final)
        else:
            location = None
        if versionControl.VersionControl().checksum_version_validate_DIAMOND(location):
            logging.info('Diamond DATABASE downloaded successfully! Consider running <checkm2 testrun> to verify everything works.')
        else:
            logging.error('Could not verify successful installation of reference database.')


    def update_database(self):
        pass

def update_checkm2():
    pass


def check_empty_dir(input_dir, overwrite=False):
    """Check the the specified directory is empty and create it if necessary."""
    if not os.path.exists(input_dir):
        make_sure_path_exists(input_dir)
    else:
        # check if directory is empty
        files = os.listdir(input_dir)
        if len(files) != 0:
            if overwrite:
                for root, dirs, files in os.walk(input_dir):
                    for f in files:
                        os.unlink(os.path.join(root, f))
                    for d in dirs:
                        shutil.rmtree(os.path.join(root, d))
            else:
                logging.error('Output directory must be empty: ' + input_dir + ' Use --force if you wish to overwrite '
                                                                               'existing directory. \n')
                sys.exit(1)

def verify_prodigal_output(prodigal_dir, ttable_dict, bin_extension):
    """Check the prodigal process was successful, matches internal list of genomes, and return list of protein files."""
    if not os.path.exists(prodigal_dir):
        # check if directory is empty
        logging.error('Error: Protein directory {} does not exist: ' + prodigal_dir + '\n')
        sys.exit(1)

    else:
        files = os.listdir(prodigal_dir)
        #check if files were generated
        if len(files) == 0:
            logging.error('Error: No protein files were generated in {}'.format(prodigal_dir))
            sys.exit(1)

        prodigal_files = []

        for f in files:
            if f.endswith('.faa'):
                protein_file = os.path.join(prodigal_dir, f)
                if os.stat(protein_file).st_size == 0:
                    logging.warning("Skipping protein file {} as it was empty.".format(protein_file))
                    del ttable_dict[f[:-4]]
                elif os.path.splitext(os.path.basename(f))[0] not in ttable_dict:
                    logging.warning("Skipping protein file {} as it was not generated by Checkm2.".format(protein_file))
                else:
                    prodigal_files.append(protein_file)
    if len(prodigal_files) == len(ttable_dict.keys()):
        return prodigal_files, ttable_dict
    else:
        logging.error('Error: List of protein files does not match internal reference.')
        sys.exit(1)

def check_if_file_exists(inputFile):
    """Check if file exists."""
    if not os.path.exists(inputFile):
        logging.error('File does not exist: ' + inputFile + '\n')
        sys.exit(1)


def check_if_dir_exists(inputDir):
    """Check if directory exists."""
    if not os.path.exists(inputDir):
        logging.error('Input directory does not exists: ' + inputDir + '\n')
        sys.exit(1)


def make_sure_path_exists(path):
    """Create directory if it does not exist."""
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logging.error('Specified path does not exist: ' + path + '\n')
            sys.exit(1)
