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
#    import fastqlib as u
#    
#    input  = u.FastQSource('/path/to/file.fastq')
#    output = u.FastQOutput('/path/to/output.fastq')
#
#    while input.next(trim_from = 5, trim_to = 75):
#        if input.entry.Q_mean > 20:
#            output.store_entry(input.entry)
#    ---------------------------------------------------------
#

import os
import sys
import stat
import gzip
import numpy
import cPickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import fastalib as u

class FastQLibError(Exception):
    pass


def compute_plot_dict_from_tiles_dict(tiles_dict, plot_dict = {'1': {}, '2': {}}):
    for pair_no in ['1', '2']:
        for tile_no in tiles_dict[pair_no]:
            if not plot_dict[pair_no].has_key(tile_no):
                plot_dict[pair_no][tile_no] = {'mean': [], 'count': [], 'std': []}
    
    for pair_no in ['1', '2']:
        for tile_no in tiles_dict[pair_no]:
            for i in range(0, 101):
                plot_dict[pair_no][tile_no]['mean'].append(numpy.mean(tiles_dict[pair_no][tile_no][i]))
                plot_dict[pair_no][tile_no]['std'].append(numpy.std(tiles_dict[pair_no][tile_no][i]))
                plot_dict[pair_no][tile_no]['count'].append(len(tiles_dict[pair_no][tile_no][i]))
                
    return plot_dict

def visualize_sequence_length_distribution(fasta_file_path, dest, title, max_seq_len = None, xtickstep = None, ytickstep = None):
    sequence_lengths = []

    fasta = u.SequenceSource(fasta_file_path)

    while fasta.next():
        if fasta.pos % 10000 == 0 or fasta.pos == 1:
            sys.stderr.write('\rReading: %s' % (big_number_pretty_print(fasta.pos)))
            sys.stderr.flush()
        sequence_lengths.append(len(fasta.seq))

    sys.stderr.write('\n')

    if not max_seq_len:
        max_seq_len = max(sequence_lengths) + (int(max(sequence_lengths) / 100.0) or 10)

    seq_len_distribution = [0] * (max_seq_len + 1)

    for l in sequence_lengths:
        seq_len_distribution[l] += 1

    fig = plt.figure(figsize = (16, 12))
    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)

    gs = gridspec.GridSpec(10, 1)

    ax1 = plt.subplot(gs[0:8])
    plt.grid(True)
    plt.subplots_adjust(left=0.05, bottom = 0.03, top = 0.95, right = 0.98)

    plt.plot(seq_len_distribution, color = 'black', alpha = 0.3)
    plt.fill_between(range(0, max_seq_len + 1), seq_len_distribution, y2 = 0, color = 'black', alpha = 0.15)
    plt.ylabel('number of sequences')
    plt.xlabel('sequence length')

    if xtickstep == None:
        xtickstep = (max_seq_len / 50) or 1

    if ytickstep == None:
        ytickstep = max(seq_len_distribution) / 20 or 1

    plt.xticks(range(xtickstep, max_seq_len + 1, xtickstep), rotation=90, size='xx-small')
    plt.yticks(range(0, max(seq_len_distribution) + 1, ytickstep), size='xx-small')
    plt.ylim(ymin = 0, ymax = max(seq_len_distribution) + (max(seq_len_distribution) / 20.0))
    plt.xlim(xmin = 0, xmax = max_seq_len)
    plt.yticks(size='xx-small')

    plt.figtext(0.5, 0.96, '%s' % (title), weight = 'black', size = 'xx-large', ha = 'center')

    ax1 = plt.subplot(gs[9])
    plt.rcParams.update({'axes.edgecolor' : 20})
    plt.grid(False)
    plt.yticks([])
    plt.xticks([])
    plt.text(0.02, 0.5, 'total: %s / mean: %.2f / std: %.2f / min: %s / max: %s'\
        % (big_number_pretty_print(len(sequence_lengths)),
           numpy.mean(sequence_lengths), numpy.std(sequence_lengths),\
           big_number_pretty_print(min(sequence_lengths)),\
           big_number_pretty_print(max(sequence_lengths))),\
        va = 'center', alpha = 0.8, size = 'x-large')


    try:
        plt.savefig(dest + '.tiff')
    except:
        plt.savefig(dest + '.png')


class Gs:
    def __init__(self, x, y):
        self.grid = gridspec.GridSpec(x, y)
        self.pointer = None
        self.current = None
    
    def next(self, p = 1):
        if self.pointer == None:
            self.pointer = 0
        else:
            self.pointer += p

        self.current = self.grid[self.pointer:self.pointer + p]
        return self.current


def visualize_qual_stats_dict(D, dest, title, split_tiles = False):
    """
    
    this how D looks like:

    D = {
            '1': {
                    'tile_1' : {'mean': [], 'std': [], 'count': []},
                    (...)
                    'tile_n' : {'mean': [], 'std': [], 'count': []},
                 },
            '2': {
                    'tile_1' : {'mean': [], 'std': [], 'count': []},
                    (...)
                    'tile_n' : {'mean': [], 'std': [], 'count': []},
                 },
        }
   
    there are two entries per pair, for each pair there are N tiles, for ever tile,
    there are mean, std and count entries, in each of those, there are 101 items,
    where nth item corresponds to the mean of all values at nth location of the
    given tile.

    dest is destination file to save the output.

    title is the title to put on the figure.
    
    """

    def get_max_count(D, p = '1', tile = None):
        if tile == None:
            return float(max([D[p][x]['count'][0] for x in D[p].keys()]))
        else:
            return float(max(D[p][tile]['count']))


    if split_tiles:
        fig = plt.figure(figsize = (30, 16))
        gs = Gs(6, 80)
    else:
        fig = plt.figure(figsize = (30, 16))
        gs = Gs(6, 16)
    
    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)
    
    
    subplots = {}
    colors = cm.get_cmap('RdYlGn', lut=256)

    tiles = ['1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108','1109', '1110', '1111','1112', '1113', '1114', '1115', '1116', 
             '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208','1209', '1210', '1211','1212', '1213', '1214', '1215', '1216', 
             '1301', '1302', '1303', '1304', '1305', '1306', '1307', '1308','1309', '1310', '1311','1312', '1313', '1314', '1315', '1316', 
             '2101', '2102', '2103', '2104', '2105', '2106', '2107', '2108','2109', '2110', '2111','2112', '2113', '2114', '2115', '2116', 
             '2201', '2202', '2203', '2204', '2205', '2206', '2207', '2208','2209', '2210', '2211','2212', '2213', '2214', '2215', '2216', 
             '2301', '2302', '2303', '2304', '2305', '2306', '2307', '2308','2309', '2310', '2311','2312', '2313', '2314', '2315', '2316']


    if split_tiles:
        m = {}
        m['1'] = get_max_count(D, '1')
        m['2'] = get_max_count(D, '2')
   
        for i in range(0, len(tiles)):
            tile = tiles[i]
           
            for _pair, _color in [('1', 'orange'), ('2', 'purple'), (None, None)]:
                if _pair:
                    subplots[tile] = {_pair: plt.subplot(gs.next(2))}
                else:
                    gs.next(1)
                    continue
                
                plt.grid(True)
                plt.subplots_adjust(left=0.02, bottom = 0.03, top = 0.95, right = 0.98)
  
                plt.xticks(range(10, 101, 10), rotation=90, size='xx-small')
                plt.ylim(ymin = 0, ymax = 42)
                plt.xlim(xmin = 0, xmax = 100)
                if _pair == '1':
                    plt.yticks(range(5, 41, 5), size='xx-small')
                else:
                    plt.yticks(range(5, 41, 5), [])

                if D[_pair].has_key(tile):
                    plt.fill_between(range(0, 101), [42 for _ in range(0, 101)], y2 = 0, color = colors(D[_pair][tile]['count'][0] / m[_pair]), alpha = 0.2)
                    subplots[tile][_pair].plot(D[_pair][tile]['mean'], color = _color, lw = 2)
                    
                    read_number_percent_dropdown = [42 * (x / get_max_count(D, _pair, tile)) for x in D[_pair][tile]['count']]
                    if not len(set(read_number_percent_dropdown)) <= 1:
                        plt.fill_between(range(0, 101), read_number_percent_dropdown, y2 = 0, color = 'black', alpha = 0.08)

                    plt.text(5, 2.5, '%s/%s :: %s' % (tile, _pair, big_number_pretty_print(int(get_max_count(D, _pair, tile)))), alpha=0.8, size = 'x-small')
                else:
                    plt.text(5, 2.5, '%s/%s :: 0' % (tile, _pair), alpha=0.8, size = 'x-small')
                    plt.fill_between(range(0, 101), [42 for _ in range(0, 101)], y2 = 0, color = colors(0 / m[_pair]), alpha = 0.2)

        Total = lambda p: big_number_pretty_print(int(sum([D[p][x]['count'][0] for x in D[p]])))

        plt.figtext(0.5, 0.97, '%s (p1: %s; p2: %s)' % (title, Total('1'), Total('2')), weight = 'black', size = 'xx-large', ha = 'center')

    else:
        m = get_max_count(D)
        
        for i in range(0, len(tiles)):
            tile = tiles[i]
            
            subplots[tile] = plt.subplot(gs.next())
            plt.grid(True)

            plt.subplots_adjust(left=0.02, bottom = 0.03, top = 0.95, right = 0.98)
  
            plt.xticks(range(10, 101, 10), rotation=90, size='xx-small')
            plt.ylim(ymin = 0, ymax = 42)
            plt.xlim(xmin = 0, xmax = 100)
            plt.yticks(range(5, 41, 5), size='xx-small')

            if D['1'].has_key(tile):
                plt.fill_between(range(0, 101), [42 for _ in range(0, 101)], y2 = 0, color = colors(D['1'][tile]['count'][0] / m), alpha = 0.2)
                subplots[tile].plot(D['1'][tile]['mean'], color = 'orange', lw = 2)
                
                read_number_percent_dropdown = [42 * (x / get_max_count(D, tile = tile)) for x in D['1'][tile]['count']]
                if not len(set(read_number_percent_dropdown)) <= 1:
                    plt.fill_between(range(0, 101), read_number_percent_dropdown, y2 = 0, color = 'black', alpha = 0.08)

                plt.text(5, 2.5, '%s :: %s' % (tile, big_number_pretty_print(int(get_max_count(D, tile = tile)))), alpha=0.5)
            else:
                plt.text(5, 2.5, '%s :: 0' % tile, alpha=0.5)
                plt.fill_between(range(0, 101), [42 for _ in range(0, 101)], y2 = 0, color = colors(0 / m), alpha = 0.2)
            if D['2'].has_key(tile):
                subplots[tile].plot(D['2'][tile]['mean'], color = 'purple', lw = 2)

        plt.figtext(0.5, 0.97, '%s (%s)' % (title, big_number_pretty_print(int(sum([D['1'][x]['count'][0] for x in D['1']])))), weight = 'black', size = 'xx-large', ha = 'center')

    try:
        plt.savefig(dest + '.tiff')
    except:
        plt.savefig(dest + '.png')


def populate_tiles_qual_dict_from_input(input_1, input_2, tiles_dict = {'1': {}, '2': {}}):
    while input_1.next():
        if input_1.p_available:
            input_1.print_percentage()
    
        if not tiles_dict['1'].has_key(input_1.entry.tile_number):
            tiles_dict['1'][input_1.entry.tile_number] = []
            for i in range(0, 101):
                tiles_dict['1'][input_1.entry.tile_number].append([])
   
        q1 = input_1.entry.process_Q_list()
        
        for i in range(0, len(q1)):
            tiles_dict['1'][input_1.entry.tile_number][i].append(q1[i])
    
    sys.stderr.write('\n') 
    
    while input_2.next():
        if input_2.p_available:
            input_2.print_percentage()

        if not tiles_dict['2'].has_key(input_2.entry.tile_number):
            tiles_dict['2'][input_2.entry.tile_number] = []
            for i in range(0, 101):
                tiles_dict['2'][input_2.entry.tile_number].append([])
        
        q2 = input_2.entry.process_Q_list()
        
        for i in range(0, len(q2)):
            tiles_dict['2'][input_2.entry.tile_number][i].append(q2[i])
    
    sys.stderr.write('\n') 
    return tiles_dict
 

def predict_file_length(file_pointer, file_path):
    file_stat = os.stat(file_path)
    file_size = file_stat[stat.ST_SIZE]

    block_length = len(''.join([file_pointer.readline() for c in range(0, 4)]))
    file_pointer.seek(0)

    if block_length:
        return file_size / block_length
    else:
        return None

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
        if key in ['__str__']:
            return None

        if key in ['trim']:
            return self.trim

        if key in ['Q_min', 'Q_mean', 'Q_std'] and self.Q_list is None:
            self.process_Q_list()

        return getattr(self, '_'.join(['process', key]))()
    
    def __init__(self, (header_line, sequence_line, optional_line, qual_scores_line), trim_from = 0, trim_to = 150, raw = False):
        self.is_valid = False

        if not header_line:
            return None

        if not header_line.startswith('@') or not optional_line.startswith('+'):
            # I'll see how stringent format checking should be. It might be a better
            # idea to raise an error than returning False.
            return None


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

    def trim(self, trim_from = 0, trim_to = 150):
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

    def next(self, trim_from = 0, trim_to = 150, raw = False):
        self.entry = FastQEntry([self.file_pointer.readline().strip() for _ in range(0, 4)], trim_from, trim_to, raw)
       
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

    def print_percentage(self, prefix = ''):
        if self.percent_read:
            sys.stderr.write('\r%s %.2d%% -- (approximate number of entries have been processed so far: %s)' % (prefix, self.percent_read, big_number_pretty_print(self.pos)))
        else:
            sys.stderr.write('\r%s (approximate number of entries have been processed so far: %s)' % (prefix, big_number_pretty_print(self.pos)))
        
        sys.stderr.flush()
        self.p_available = False
    
    def close(self):
        self.file_pointer.close()

