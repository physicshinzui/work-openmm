#!/bin/bash 
set -eu
INITIAL_STRUCUTURE="../../sample_pdb/capped_tri_ala.pdb"
#INITIAL_STRUCUTURE="pdb_name_is_here"

python ../../../src/runner.py $INITIAL_STRUCUTURE
