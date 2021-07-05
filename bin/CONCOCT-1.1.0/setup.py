#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
from sys import platform
from distutils.core import Extension
import numpy as np

try:
    from Cython.Distutils import build_ext
except ImportError:
    print("You need to have Cython installed on your system to run setup.py. Sorry!")
    sys.exit()

version = '1.1.0'

include_dirs_for_concoct = [np.get_include(), '/opt/local/include/']     

extra_compile_args = ['-O3','-std=c99']
# System clang on MacOS does not recognize the -fopenmp argument
if platform != 'darwin':
    extra_compile_args = ['-fopenmp'] + extra_compile_args

setup(name='concoct',
      version=version,
      description="Clustering cONtigs with COverage and ComposiTion",
      long_description="""Concoct is a program that combines three types of
      information - sequence composition, coverage across multiple sample,
      and read-pair linkage - to automatically bin metagenomic contigs
      into genomes. """,
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg, Christopher Quince, Anders Andersson, Ino de Bruijn',
      author_email='binni@binnisb.com',
      maintainer='Johannes Alneberg',
      maintainer_email='johannes.alneberg@scilifelab.se',
      url='https://github.com/BinPro/CONCOCT',
      license='FreeBSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["bin/concoct","bin/concoct_refine", "scripts/cut_up_fasta.py", "scripts/concoct_coverage_table.py",
               "scripts/merge_cutup_clustering.py", "scripts/extract_fasta_bins.py"],
      include_package_data=True,
      zip_safe=False,
      cmdclass = {'build_ext': build_ext},
      ext_modules = [
                    Extension("vbgmm", sources=["./c-concoct/vbgmm.pyx", "./c-concoct/c_vbgmm_fit.c"],
                                libraries =['gsl', 'gslcblas', 'gomp'], include_dirs=include_dirs_for_concoct,
                                extra_compile_args = extra_compile_args)
                    ],
      install_requires=['cython>=0.19.1',
                        'numpy>=1.7.1',
                        'scipy>=0.12.0',
                        'pandas>=0.11.0',
                        'biopython>=1.62b',
                        'scikit-learn>=0.13.1',
                        'nose'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
