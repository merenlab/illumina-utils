# -*- coding: utf-8 -*-
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
#    import fastq as u
#    
#    input  = u.FastQSource('/path/to/file.fastq')
#    output = u.FastQOutput('/path/to/output.fastq')
#
#    while input.next(trim_to = 75):
#        if input.entry.Q_mean > 20:
#            output.store(input.entry)
#    ---------------------------------------------------------
#

import os
import sys
import stat
import numpy

def predict_file_length(file_pointer, file_path):
    file_stat = os.stat(file_path)
    file_size = file_stat[stat.ST_SIZE]

    block_length = len(''.join([file_pointer.readline() for c in range(0, 4)]))
    file_pointer.seek(0)

    return file_size / block_length

def big_number_pretty_print(n):
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()

    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

class FastQEntry:
    def __getattr__(self, key):
        # on demand quality processing.
        if key in ['Q_min', 'Q_mean', 'Q_std'] and self.Q_list is None:
            self.process_Q_list()

        return getattr(self, '_'.join(['process', key]))()
    
    def __init__(self, (header_line, sequence_line, optional_line, qual_scores_line), trim_to = 150):
        self.is_valid = False

        if not header_line:
            return None

        if not header_line.startswith('@') or not optional_line.startswith('+'):
            # I'll see how stringent format checking should be. It might be a better
            # idea to raise an error than returning False.
            return None

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

        self.sequence    = sequence_line[0:self.trim_to]
        self.qual_scores = qual_scores_line[0:self.trim_to]
        self.optional    = optional_line[1:]

        self.Q_list = None

        self.is_valid = True

    def process_Q_list(self):
        self.Q_list = [ord(q) - 33 for q in self.qual_scores]
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

class FastQOutput:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_pointer = open(file_path, 'w')

    def store(self, e):
        header_line = ':'.join([e.machine_name, e.run_id, e.flowcell_id, 
                                e.lane_number, e.tile_number, e.x_coord,
                                e.y_coord + ' ' + e.pair_no if e.pair_no else e.y_coord,
                                'Y' if e.quality_passed else 'N', e.control_bits_on, e.index_sequence])

        self.file_pointer.write('@' + header_line + '\n')
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
    def __init__(self, file_path):
        self.pos = 0

        self.file_pointer = open(file_path)
        self.file_length = predict_file_length(self.file_pointer, file_path)

        self.percent_step = self.file_length / 1000
        self.percent_read = None
        self.p_available = False
        self.percent_counter = 0

    def next(self, trim_to = 150):
        self.entry = FastQEntry([self.file_pointer.readline().strip() for _ in range(0, 4)], trim_to)
       
        if not self.entry.is_valid:
            return False

        self.pos += 1

        if self.pos == 1 or self.pos % 1000 == 0:
            self.p_available = True
            self.percent_read = self.pos * 100 / self.file_length

        return True

    def reset(self):
        self.pos = 0
        self.file_pointer.seek(0)
    
    def close(self):
        self.file_pointer.close()

    def print_percentage(self):
        sys.stderr.write('\r%.2d%% -- (approximate number of entries have been processed so far: %s)' % (self.percent_read, big_number_pretty_print(self.pos)))
        sys.stderr.flush()
        self.p_available = False

