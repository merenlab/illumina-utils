# -*- coding: utf-8 -*-
#
# Copyright (C) 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#
#
#
# Evolving Python utility library for Illumina stuff (for CASAVA 1.7+ pipeline).
#
# Questions, suggestions: A. Murat Eren <a.murat.eren / gmail>
#
# fastq.py -- this file is supposed to be a common library to read/write fastq entries and
# process qual scores that might be required for various operations.
#
#
# Example usages of this library:
#
#    (trim every read to 75 and then filter the ones that have
#     mean PHRED qual score above Q20 into another file).
#    ---------------------------------------------------------
#    import fastqlib as u
#    
#    input  = u.FastQSource('/path/to/file.fastq')
#    output = u.FastQOutput('/path/to/output.fastq')
#
#    while input.next(trim_from = 5, trim_to = 75):
#        if input.entry.Q_mean > 20:
#            output.store_entry(input.entry)
#    ---------------------------------------------------------


import os
import sys
import gzip
import numpy

from IlluminaUtils.utils.helperfunctions import big_number_pretty_print
from IlluminaUtils.utils.helperfunctions import predict_file_length


class FastQLibError(Exception):
    pass


class FastQEntry:
    def __getattr__(self, key):
        # on demand quality processing.
        if key in ['__str__']:
            return None

        if key in ['trim']:
            return self.trim

        if key in ['Q_min', 'Q_mean', 'Q_std'] and self.Q_list is None:
            self.process_Q_list()

        return getattr(self, '_'.join(['process', key]))()
    
    def __init__(self, (header_line, sequence_line, optional_line, qual_scores_line), trim_from = 0, trim_to = sys.maxint, raw = False, CASAVA_version = 1.8):
        self.is_valid = False

        if not header_line:
            return None

        if not header_line.startswith('@') or not optional_line.startswith('+'):
            # I'll see how stringent format checking should be. It might be a better
            # idea to raise an error than returning False.
            return None

        self.CASAVA_version = CASAVA_version

        self.header_line = header_line[1:]

        if not raw:
            header_fields = header_line[1:].split(':')

            self.machine_name   = header_fields[0]
            self.run_id         = header_fields[1]
            self.flowcell_id    = header_fields[2]
            self.lane_number    = header_fields[3]
            self.tile_number    = header_fields[4]
            self.x_coord        = header_fields[5]
            if len(header_fields[6].split()):
                self.y_coord    = header_fields[6].split()[0]
                self.pair_no    = header_fields[6].split()[1]
            else:
                self.y_coord    = header_fields[6]
                self.pair_no    = None
            self.quality_passed = header_fields[7] == 'Y'
            self.control_bits_on = header_fields[8]
            self.index_sequence  = header_fields[9]

        self.trim_to = trim_to
        self.trim_from = trim_from

        self.sequence    = sequence_line[self.trim_from:self.trim_to]
        self.qual_scores = qual_scores_line[self.trim_from:self.trim_to]
        self.optional    = optional_line[1:]

        self.Q_list = None

        self.is_valid = True

    def process_Q_list(self):
        if self.CASAVA_version >= 1.8:
            self.Q_list = [ord(q) - 33 for q in self.qual_scores]
        if self.CASAVA_version < 1.8:
            self.Q_list = [ord(q) - 64 for q in self.qual_scores]
 
        return self.Q_list
    
    def process_Q_min(self):
        self.Q_min = numpy.min(self.Q_list)
        return self.Q_min
        
    def process_Q_mean(self):
        self.Q_mean = numpy.mean(self.Q_list)
        return self.Q_mean

    def process_Q_std(self):
        self.Q_std  = numpy.std(self.Q_list)
        return self.Q_std

    def trim(self, trim_from = 0, trim_to = sys.maxint):
        self.trim_to = trim_to
        self.trim_from = trim_from
        self.sequence    = self.sequence[trim_from:self.trim_to]
        self.qual_scores = self.qual_scores[trim_from:self.trim_to]


class FileOutput(object):
    def __init__(self, file_path, compressed = False):
        self.file_path = file_path
        
        self.compressed = compressed
        
        if self.compressed:
            self.file_pointer = gzip.open(file_path, 'w')
        else:
            self.file_pointer = open(file_path, 'w')

    def write(self, c):
        self.file_pointer.write(c)

    def close(self):
        self.file_pointer.close()


class FastQOutput(FileOutput):
    def __init__(self, file_path, compressed = False):
        # lets don't overwrite anything and protect users from their
        # mistakes. Nope; I didn't overwrite my millions of sequences.
        # Not at all.
        if os.path.exists(file_path):
            raise FastQLibError, 'Output file exists: "%s"' % file_path

        super(FastQOutput, self).__init__(file_path, compressed)

    def store_entry(self, e):
        self.file_pointer.write('@' + e.header_line + '\n')
        self.file_pointer.write(e.sequence + '\n')
        self.file_pointer.write('+' + e.optional + '\n')
        self.file_pointer.write(e.qual_scores + '\n')


class FastQSource:
    """
    Class to iterate entries in a fastq file.

    ---------------------------------------------------------
    import fastq as u
    
    input  = u.FastQSource('/path/to/file.fastq')

    while input.next(trim_to = 75):
        print input.entry.Q_max, input.entry.sequence
    ---------------------------------------------------------
    """
    def __init__(self, file_path, compressed = False):
        if not os.path.exists(file_path):
            raise FastQLibError, 'Missing input file: "%s"' % file_path
        
        self.pos = 0

        self.compressed = compressed
        if self.compressed:
            self.file_pointer = gzip.open(file_path)
            self.file_length = None
        else:
            self.file_pointer = open(file_path)
            self.file_length = predict_file_length(self.file_pointer, file_path)

        self.percent_step = 1000
        self.percent_read = None
        self.p_available = False
        self.percent_counter = 0
        
        # FIXME: check for actual CASAVA version
        self.CASAVA_version = 1.8

    def next(self, trim_from = 0, trim_to = sys.maxint, raw = False):
        self.entry = FastQEntry([self.file_pointer.readline().strip() for _ in range(0, 4)], trim_from, trim_to, raw, CASAVA_version = self.CASAVA_version)
       
        if not self.entry.is_valid:
            return False

        self.pos += 1

        if self.pos == 1 or self.pos % 1000 == 0:
            self.p_available = True
            if self.file_length:
                self.percent_read = self.pos * 100 / self.file_length
            else:
                self.percent_read = None

        return True

    def reset(self):
        self.pos = 0
        self.file_pointer.seek(0)
    
    def close(self):
        self.file_pointer.close()

    def print_percentage(self, prefix = '', postfix = ''):
        if self.percent_read:
            sys.stderr.write('\r%s %.2d%% -- (num pairs processed: %s) %s' % (prefix, self.percent_read, big_number_pretty_print(self.pos), postfix))
        else:
            sys.stderr.write('\r%s (num pairs processed: %s) %s' % (prefix, big_number_pretty_print(self.pos), postfix))
        
        sys.stderr.flush()
        self.p_available = False


if __name__ == "__main__":
    pass
