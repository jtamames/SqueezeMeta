from checkm2.defaultValues import DefaultValues
from checkm2 import fileManager
from checkm2 import version

import os
import hashlib
from packaging import version as v_compare
import pandas as pd
import json
import sys
import logging


class VersionControl():

    def __init__(self):
        self.version = version.__version__

    def checksum_test_genomes(self):

        ''' Hardcoded until more robust testing is set up'''

        checksums = {0: 'f097f8a99ad692f26b758387c5374f426017eaf9251a41f35d5bf77e1bd8bc5a',
                     1: 'a91663bb71937d539e3f55cedd0709be1cd4aaa3e1c95b30b985691108977acc',
                     2: '53d80449105af334d523cda16803fc1ae9247f3f463468eed857f6ca466bfee2'}

        TEST1 = os.path.join(DefaultValues.TESTRUN_GENOMES, 'TEST1.tst')
        TEST2 = os.path.join(DefaultValues.TESTRUN_GENOMES, 'TEST2.tst')
        TEST3 = os.path.join(DefaultValues.TESTRUN_GENOMES, 'TEST3.tst')

        checksum_results = []

        for idx, tg in enumerate([TEST1, TEST2, TEST3]):
            readable_hash = self.__calculate_checksum(tg)
            checksum_results.append(readable_hash == checksums[idx])

        if False not in checksum_results:
            return True
        else:
            return False


    def __calculate_checksum(self, query_file, chunk=False):
        if not chunk:
            with open(query_file, "rb") as f:
                bytes = f.read()  # read entire file as bytes
                return hashlib.sha256(bytes).hexdigest()
        else:
            #for larger files
            sha256_hash = hashlib.sha256()
            with open(query_file, "rb") as f:
                # Read and update hash string value in blocks of 4K
                for byte_block in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(byte_block)
                return sha256_hash.hexdigest()

    def __validateVersion(self, query, ref):
        return v_compare.parse(str(query)) >= v_compare.parse(str(ref))

    def return_highest_compatible_DB_version(self):
        version_hashes = os.path.join(DefaultValues.VERSION_PATH, 'version_hashes_{}.json'.format(version.__version__))
        try:
            version_hashes = pd.read_json(version_hashes)
        except Exception as e:
            logging.error('Could not open correct hash version file: {}'.format(e))
            sys.exit(1)
        DB_version = version_hashes[version_hashes['type'] == 'DIAMONDDB'].sort_values(by='version', ascending=False)
        DB_version['Valid'] = DB_version.apply(
            lambda row: self.__validateVersion(row['version'], row['incompatible_below_checkm2ver']), axis=1)

        return DB_version[DB_version['Valid'] == True]['version'].values[0], DB_version[DB_version['Valid'] == True]['DOI'].values[0]


    def checksum_version_validate(self):
        '''Runs each time to ensure all models, definitions and pickled files are congruent with current CheckM2 version'''

        version_hashes = os.path.join(DefaultValues.VERSION_PATH, 'version_hashes_{}.json'.format(version.__version__))
        try:
            version_hashes = pd.read_json(version_hashes)
        except Exception as e:
            logging.error('Could not open correct hash version file: {}'.format(e))
            sys.exit(1)


        for file in DefaultValues.EXTERNAL_FILES_TO_VERIFY:
            filehash = self.__calculate_checksum(file)
            if filehash in version_hashes['sha256'].values:
                cutoff_version = version_hashes[version_hashes["sha256"] == filehash]['incompatible_below_checkm2ver']
                if self.__validateVersion(self.version, cutoff_version):
                    return True
                else:
                    return False
            else:
                logging.error('Fault was with file: {}'.format(file))
                logging.error('Hash was {}'.format(filehash))
                logging.error('One of the files CheckM2 relies on has an incorrect checksum. Please re-download CheckM2.')
                sys.exit(1)

    def checksum_version_validate_DIAMOND(self, location=None):
        '''Runs to ensure DIAMOND database has correct checksum and is congruent with current CheckM2 version'''

        version_hashes = os.path.join(DefaultValues.VERSION_PATH,
                                      'version_hashes_{}.json'.format(version.__version__))
        try:
            version_hashes = pd.read_json(version_hashes)
        except Exception as e:
            logging.error('Could not open correct hash version file: {}'.format(e))
            sys.exit(1)

        if not location:
            database_dir = fileManager.DiamondDB().get_DB_location()
        else:
            database_dir = location

        dbhash = self.__calculate_checksum(database_dir, True)
        

        if dbhash in version_hashes['sha256'].values:
            cutoff_version = version_hashes[version_hashes["sha256"] == dbhash]['incompatible_below_checkm2ver'].values[-1]
            return self.__validateVersion(self.version, cutoff_version)
        else:
            logging.error('One of the files CheckM2 relies on has an incorrect checksum. Please re-download CheckM2.')
            sys.exit(1)


