#!/bin/bash
# this file is to install Short-Pair in the current folder.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
chmod 755 Short-Pair.py
chmod 755 get_hmm.sh
chmod 755 hmmer3_pipeline_missing_end.sh
