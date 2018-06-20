#!/bin/bash
source info.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd demultiplexing
rm -rf output
mkdir output

INFO "Interleaving FASTQ files ..."
iu-interleave-fastq -1 r1.fastq -2 r2.fastq -o output/interleaved.fastq

INFO "De-interleaving FASTQ files ..."
iu-deinterleave-fastq output/interleaved.fastq -1 output/deinterleaved_R1.fastq -2 output/deinterleaved_R2.fastq

INFO "Done."
