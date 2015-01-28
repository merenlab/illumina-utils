#!/bin/bash
source info.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd v6-complete-overlap
rm -rf output
mkdir output

INFO "Generating the config file from sample descripton"
iu-gen-configs sample-description.txt -o output --r1-prefix 'A[AC][CG]CG[AC].[AG]AACCT[CT]A.C' --r2-prefix '^....C[AG][AG]CCATGCA.CACCT'

INFO "Running complete overlap analysis for V6 ..."
iu-merge-pairs output/E.coli.ini --marker-gene-stringent --retain-only-overlap

INFO "Retaining only reads with 0 mismatches"
iu-filter-merged-reads output/E.coli_MERGED --max-mismatches 0


INFO "Done."
