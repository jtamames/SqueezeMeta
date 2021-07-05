# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 2013

@author: Johannes Alneberg
"""
from __future__ import print_function
import os
import sys
import logging

import pandas as p
import numpy as np


class Output(object):
    """
    Class to print out result information to their files
    """
    CONCOCT_PATH = None
    ARGS_FILE = None
    PCA_FILE_BASE = None
    PCA_COMPONENTS_FILE_BASE = None
    FLOAT_FORMAT = '%1.8e'
    INT_FORMAT = '%d'

    @classmethod
    def __init__(self,basename,args):
        """
        Set output params and create output folders and bic.csv and args.txt
        """
        if os.path.isdir(basename):
            if basename[-1] == '/':
                self.CONCOCT_PATH = basename
            else:
                self.CONCOCT_PATH = basename + '/'
        elif basename[-1] == '/':
            basename_path = os.path.abspath(basename)
            os.mkdir(basename_path)
            self.CONCOCT_PATH = basename_path +'/'
        else:
            basename_path = os.path.abspath(basename)
            self.CONCOCT_PATH = basename+'_'

        self.ARGS_FILE = self.CONCOCT_PATH + "args.txt"
        self.ORIGINAL_FILE_BASE = self.CONCOCT_PATH + "original_data_gt{0}.csv"
        self.PCA_FILE_BASE = self.CONCOCT_PATH + \
            "PCA_transformed_data_gt{0}.csv"
        self.ASSIGN_FILE_BASE = self.CONCOCT_PATH + \
            "clustering_gt{0}.csv"
        self.PCA_COMPONENTS_FILE_BASE = self.CONCOCT_PATH + \
            "PCA_components_data_gt{0}.csv"
        self.LOG_FILE_BASE = self.CONCOCT_PATH + 'log.txt'

        # Reset any previous logging handlers, see:
        # https://stackoverflow.com/a/49202811
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(
            filename=self.LOG_FILE_BASE,
            level=logging.INFO,
            filemode='w', # Overwrites old log file
            format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
            )

        logging.info("Results created at {0}".format(
            os.path.abspath(self.CONCOCT_PATH)))


        print("Up and running. Check {0} for progress".format(
                        os.path.abspath(self.LOG_FILE_BASE)
                    ), file=sys.stderr)

        #Write header to bic.csv
        with open(self.ARGS_FILE,"w+") as fh:
            print(args, file=fh)

    @classmethod
    def write_pca(self, transform, threshold, index):
        transform_df = p.DataFrame(transform, index=index)
        transform_df.to_csv(
            self.PCA_FILE_BASE.format(threshold),
            float_format=self.FLOAT_FORMAT,
            index_label="contig_id"
            )
        logging.info('Wrote PCA transformed file.')

    @classmethod
    def write_assign(self, assign, threshold, index):
        transform_df = p.DataFrame(assign, index=index)
        transform_df.columns=['cluster_id']
        transform_df.to_csv(
            self.ASSIGN_FILE_BASE.format(threshold),
            index_label="contig_id"
            )
        logging.info('Wrote assign file.')


    @classmethod
    def write_pca_components(self, components, threshold):
        np.savetxt(
            self.PCA_COMPONENTS_FILE_BASE.format(threshold),
            components,
            fmt=self.FLOAT_FORMAT,
            delimiter=","
        )
        logging.info('Wrote PCA components file.')

    @classmethod
    def write_original_data(self,original,threshold):
        original.to_csv(self.ORIGINAL_FILE_BASE.format(threshold), float_format=self.FLOAT_FORMAT)
        logging.info('Wrote original filtered data file.')
