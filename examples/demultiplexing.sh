#!/bin/bash
source info.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd demultiplexing
rm -rf output
mkdir output

INFO "Running demultiplexing ..."
iu-demultiplex -s barcode_to_sample.txt --r1 r1.fastq --r2 r2.fastq --index index.fastq -o output/

INFO "Generating config files for each demultiplexed sample"

iu-gen-configs output/00_DEMULTIPLEXING_REPORT

INFO "Done."
