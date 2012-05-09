# -*- coding: utf-8 -*-

import os
import sys
import cPickle
import hashlib
import tempfile
import operator
import subprocess
import ConfigParser

try:
    import Levenshtein as l
except:
    print '''
    You need Levenshtein module installed to run this software.

    Here is a fast implementation of Levenshtein distance for Python:

        http://code.google.com/p/pylevenshtein/

'''
    sys.exit(-1)

import fastqlib as u
import fastalib as f

from runconfiguration import RunConfiguration
from runconfiguration import ConfigError

conv_dict = {'A': 'T',
             'T': 'A',
             'C': 'G',
             'G': 'C',
             'N': 'N'}

NumberOfConflicts = lambda s: sum([True for n in s if n in conv_dict])
ConflictPositions = lambda s: [i for i in range(0, len(s)) if s[i] in conv_dict]

def merge_two_sequences(seq1, seq2, min_overlap_size):
    smallest, ind = sys.maxint, 0
    for i in range(min_overlap_size, len(seq1)):
        d = l.distance(seq1[-i:], seq2[:i])
        if d <= smallest:
            smallest = d
            ind = i

    beg = seq1[0:len(seq1) - ind].lower()
    overlap_1 = seq1[len(seq1) - ind:].lower()
    overlap_2 = seq2[:ind].lower()
    end = seq2[ind:].lower()
    overlap = ''
    
    mismatches = 0
    for i in range(0, len(overlap_1)):
        if overlap_1[i] != overlap_2[i]:
            mismatches += 1
            if reverse:
                overlap += overlap_2[i].upper()
            else:
                overlap += overlap_1[i].upper()
        else:
            overlap += overlap_1[i]

    return (beg, overlap, end, mismatches)


def reverse_complement(seq):
    return ''.join(reversed([conv_dict[n] for n in seq]))

def reverse(seq):
    return ''.join(reversed(seq))

def complement(seq):
    return ''.join([conv_dict[n] for n in seq])


def main(config, output_file_prefix, skip_qual_dicts = False, min_overlap_size = 15):
   
    quality_dicts = {-1: {'1': {}, '2': {}}, 0: {'1': {}, '2': {}}, 1: {'1': {}, '2': {}}, 2: {'1': {}, '2': {}}, 3: {'1': {}, '2': {}}, 4: {'1': {}, '2': {}}, 5: {'1': {}, '2': {}}}
            
    output = f.FastaOutput(os.path.join(config.output_directory, output_file_prefix + '_MERGED' ))
    failed = f.FastaOutput(os.path.join(config.output_directory, output_file_prefix + '_FAILED' ))
    
    ########################################################################################################################
    # merging..
    ########################################################################################################################

    for index in range(0, len(config.lane_1)):
        try:
            F = lambda x: os.path.basename(x).split('.')[0]
            input_1 = u.FastQSource(config.lane_1[index], compressed = False)
            input_2 = u.FastQSource(config.lane_2[index], compressed = False)

        except u.FastQLibError, e:
            print "FastQLib is not happy.\n\n\t", e, "\n"
            sys.exit()

        while input_1.next() and input_2.next():
            if input_1.p_available: 
                input_1.print_percentage('[Merging %d of %d]' % (index + 1, len(config.lane_1)))

            beg, overlap, end, mismatches = merge_two_sequences(input_1.entry.sequence, reverse_complement(input_2.entry.sequence)[:-1], min_overlap_size)

            merged_sequence = '%s%s%s' % (beg, overlap, end)

            p = 1.0 * mismatches / len(overlap)

            if p > 0.3:
                failed.write_id('%s|o/m:%f|mismatches:%d' % (input_1.entry.header_line, p, mismatches))
                failed.write_seq(merged_sequence, split = False)
            else:
                output.write_id('%s|o/m:%f|mismatches:%d' % (input_1.entry.header_line, p, mismatches))
                output.write_seq(merged_sequence, split = False)

            if not skip_qual_dicts:
                ################ quality dicts associated stuff ####################
                
                minimum_length = (len(input_1.entry.sequence) + len(input_2.entry.sequence)) - min_overlap_size
                
                if len(merged_sequence) > minimum_length or p > 0.3:
                    ind = -1
                elif mismatches >= 5:
                    ind = 5
                else:
                    ind = mismatches

                tiles_dict = quality_dicts[ind]

                if not tiles_dict['1'].has_key(input_1.entry.tile_number):
                    tiles_dict['1'][input_1.entry.tile_number] = []
                    for i in range(0, 101):
                        tiles_dict['1'][input_1.entry.tile_number].append([])
   
                q1 = input_1.entry.process_Q_list()
                
                for i in range(0, len(q1)):
                    tiles_dict['1'][input_1.entry.tile_number][i].append(q1[i])

                if not tiles_dict['2'].has_key(input_2.entry.tile_number):
                    tiles_dict['2'][input_2.entry.tile_number] = []
                    for i in range(0, 101):
                        tiles_dict['2'][input_2.entry.tile_number].append([])
                
                q2 = input_2.entry.process_Q_list()
                
                for i in range(0, len(q2)):
                    tiles_dict['2'][input_2.entry.tile_number][i].append(q2[i])
                
                ################ / quality dicts associated stuff ####################
 
        print
        input_1.close()
        input_2.close()

    output.close()
    failed.close()

    if not skip_qual_dicts:
        ################ quality dicts associated stuff ####################
        for ind in quality_dicts:
            qual_dict = quality_dicts[ind]
            if ind == -1:
                cPickle.dump(qual_dict, open(os.path.join(config.output_directory, output_file_prefix + '_MERGED_BAD_LENGTH_QUAL_DICT'), 'w'))
            else:
                cPickle.dump(qual_dict, open(os.path.join(config.output_directory, output_file_prefix + '_MERGED_%d_MISMATCH_QUAL_DICT' % ind), 'w'))
        ################ / quality dicts associated stuff ####################


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Merge Overlapping Paired-End Illumina Reads')
    parser.add_argument('user_config', metavar = 'CONFIG_FILE',
                                        help = 'User configuration to run')
    parser.add_argument('output_file_prefix', metavar = 'OUTPUT_FILE_PREFIX',
                                        help = 'Output file prefix (which will be used as a prefix\
                                                for files that appear in output directory)')
    parser.add_argument('--min-overlap-size', metavar = 'MIN_OVERLAP_SIZE', default = 15, type = int,
                                        help = 'Minimum number of nucleotides expected to overlap (default: 15)')
    parser.add_argument('--skip-qual-dicts', action = 'store_true', default = False,
                                        help = 'When set, qual ditcs will not be computed.')

    args = parser.parse_args()
    user_config = ConfigParser.ConfigParser()
    user_config.read(args.user_config)

    try: 
        config = RunConfiguration(user_config)
    except ConfigError, e:
        print "There is something wrong with the config file. This is what we know: \n\n", e
        print
        sys.exit()

    sys.exit(main(config, args.output_file_prefix, args.skip_qual_dicts, args.min_overlap_size))

