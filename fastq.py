# -*- coding: utf-8 -*-
#
# Evolving Python utility library for Illumina stuff (for CASAVA 1.7+ pipeline).
#
# Questions, suggestions: A. Murat Eren <a.murat.eren / gmail>
#
# fastq.py -- this file is supposed to be a common library to read/write fastq entries and
# process qual scores that might be required for various operations.
#

import os
import sys
import stat
import numpy

def predict_file_length(f_path):
    f_stat = os.stat(f_path)
    f_size = f_stat[stat.ST_SIZE]

    f = open(f_path)
    line_length = len(''.join([f.readline() for c in range(0, 4)]))
    f.close()

    return f_size / line_length

def big_number_pretty_print(n):
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    if ret[0] == ',':
        return ''.join(ret[1:])
    else:
        return ''.join(ret)

class SequenceSource:
    """
    Class to iterate entries in a fastq file. An example usage:

    ---------------------------------------------------------
    import fastq as u
    
    fastq_entry = u.SequenceSource('/path/to/file.fastq', trim_to = 80)

    while fastq_entry.next(process_quals = True):
        if fastq_entry.Q_mean > 20:
            # do something with fastq_entry
            print fastq_entry.sequence       # <-- only if you are very bored
    ---------------------------------------------------------
    """
    def __init__(self, f_path, trim_to = 100):
        self.pos = 0
        self.trim_to = trim_to

        self.file_pointer = open(f_path)
        self.file_length = predict_file_length(f_path)

        self.percent_step = self.file_length / 1000
        self.percent_read = None
        self.p_available = False
        self.percent_counter = 0

    def next(self, process_quals = False):
        header, sequence, optional, qual_scores = [self.file_pointer.readline().strip() for _ in range(0, 4)]
       
        if not header:
            return False

        if not header.startswith('@') or not optional.startswith('+'):
            # I'll see how stringent format checking should be. It might be a better
            # idea to raise an error than returning False.
            return False

        fields = header[1:].split(':')

        self.machine_name   = fields[0]
        self.run_id         = fields[1]
        self.flowcell_id    = fields[2]
        self.lane_number    = fields[3]
        self.tile_number    = fields[4]
        self.x_coord        = fields[5]
        if len(fields[6].split()):
            self.y_coord    = fields[6].split()[0]
            self.pair_no    = fields[6].split()[1]
        else:
            self.y_coord    = fields[6]
            self.pair_no    = None
        self.quality_passed = fields[7] == 'Y'
        self.control_bits_on = fields[8]
        self.index_sequence  = fields[9]

        self.sequence    = sequence[0:self.trim_to]
        self.qual_scores = qual_scores[0:self.trim_to]
        self.optional    = optional[1:]

        self.Q_list = [ord(q) - 33 for q in self.qual_scores] if process_quals else None
        # FIXME: not all of these are always necessary. some smarts is required here.
        self.Q_min  = numpy.min(self.Q_list) if process_quals else None
        self.Q_mean = numpy.mean(self.Q_list) if process_quals else None
        self.Q_std  = numpy.std(self.Q_list) if process_quals else None

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

