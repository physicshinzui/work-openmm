#!/bin/bash 
set -eu
INITIAL_STRUCUTURE="chignolin.pdb"
#INITIAL_STRUCUTURE="pdb_name_is_here"

python ../../../src/nvt_in_water.py $INITIAL_STRUCUTURE
