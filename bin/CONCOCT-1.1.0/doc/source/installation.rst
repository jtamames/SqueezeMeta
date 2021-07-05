Installation
============

With Bioconda [Recommended]
---------------------------

The easiest and recommended way to install concoct is through `Bioconda <https://bioconda.github.io/>`_ and `conda <https://docs.conda.io/en/latest/>`_ in an isolated environment:

::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    conda create -n concoct_env python=3 concoct


Note for Mac OSX users
~~~~~~~~~~~~~~~~~~~~~~

Currently concoct on Mac OSX can only run in single threaded mode, which drastically increases the runtime.
However, the Mac OSX installation of concoct can still be useful for testing purposes and is possible to install through conda as shown above.


Manual Installation
-------------------

The conda installation should be enough for most users.
However, if you want to modify the source code, a manual installation might be needed.
An example of a manual installation on an Ubuntu system can be seen in the `Travis CI config file <https://github.com/BinPro/CONCOCT/blob/develop/.travis.yml>`_.



Using Docker
------------

We provide a Docker image:
binpro/concoct\_latest which contains CONCOCT and its dependencies for a basic workflow.

Assuming DOcker is installed, the following command will then download the image from the Docker image
index, map the Data folder to the image and log you into the docker image.

::

    docker run -v /home/USER/Data:/opt/Data -i -t binpro/concoct_latest bash

To test concoct you can then do:

::

    $ cd /opt/CONCOCT_latest
    $ nosetests

Which should execute all tests without errors.


Other Programs Needed
---------------------

-  For assembly, use your favourite. Here is a good one:

   -  `Megahit <https://github.com/voutcn/megahit>`__


-  To create the input bam files, choose your favourite aligner, for example bowtie2 or bwa.
-  For validation of clustering using single-copy core genes we recommend using:

   -  `CheckM <https://github.com/Ecogenomics/CheckM>`__
