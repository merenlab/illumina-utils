#!/bin/bash
source info.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd v6-complete-overlap
rm -rf output
mkdir output

INFO "Generating the config file from sample descripton"
iu-gen-configs sample-description.txt -o output 

INFO "Running complete overlap analysis with 0 mismatches ..."
iu-merge-pairs output/E.coli.ini --marker-gene-stringent --retain-only-overlap --max-num-mismatches 0

INFO "Trimming V6 primers from merged reads ..."
iu-trim-V6-primers output/E.coli_MERGED

INFO "Done. Please see v6-complete-overlap/output/E.coli_MERGED_V6_PRIMERS_REMOVED"
