# -*- coding: utf-8 -*-
# v.120711

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the docs/COPYING file.

import os
import sys
import numpy


class FastaOutput:
    def __init__(self, output_file_path):
        self.output_file_path = output_file_path
        self.output_file_obj = open(output_file_path, 'w')

    def store(self, entry, split = True):
        self.write_id(entry.id)
        self.write_seq(entry.seq)

    def write_id(self, id):
        self.output_file_obj.write('>%s\n' % id)

    def write_seq(self, seq, split = True):
        self.output_file_obj.write('%s\n' % self.split(seq) if split else seq)

    def split(self, sequence, piece_length = 80):
        ticks = range(0, len(sequence), piece_length) + [len(sequence)]
        return '\n'.join([sequence[ticks[x]:ticks[x + 1]] for x in range(0, len(ticks) - 1)])

    def close(self):
        self.output_file_obj.close()

class SequenceSource:
    def __init__(self, fasta_file_path, lazy_init = True):
        self.fasta_file_path = fasta_file_path
        self.pos = 0
        self.id  = None
        self.seq = None
        self.file_pointer = open(self.fasta_file_path)
        self.file_pointer.seek(0)
        
        if lazy_init:
            self.total_seq = None
        else:
            self.total_seq = len([l for l in self.file_pointer.readlines() if l.startswith('>')])

    def next(self):
        self.id = self.file_pointer.readline()[1:].strip()
        self.seq = None

        sequence = ''
        while 1:
            line = self.file_pointer.readline()
            if not line:
                if len(sequence):
                    self.seq = sequence
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            sequence += line.strip()

        self.seq = sequence
        self.pos += 1
        return True

    def reset(self):
        self.pos = 0
        self.id  = None
        self.seq = None
        self.file_pointer.seek(0)

    def visualize_sequence_length_distribution(self, title, dest = None, max_seq_len = None, xtickstep = None, ytickstep = None):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        sequence_lengths = []
    
        self.reset()
    
        while self.next():
            if self.pos % 10000 == 0 or self.pos == 1:
                sys.stderr.write('\r[fastalib] Reading: %s' % (self.pos))
                sys.stderr.flush()
            sequence_lengths.append(len(self.seq))
        
        self.reset()
    
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
        plt.yticks(range(0, max(seq_len_distribution) + 1, ytickstep),
                   [y for y in range(0, max(seq_len_distribution) + 1, ytickstep)],
                   size='xx-small')
        plt.xlim(xmin = 0, xmax = max_seq_len)
        plt.ylim(ymin = 0, ymax = max(seq_len_distribution) + (max(seq_len_distribution) / 20.0))
    
        plt.figtext(0.5, 0.96, '%s' % (title), weight = 'black', size = 'xx-large', ha = 'center')
    
        ax1 = plt.subplot(gs[9])
        plt.rcParams.update({'axes.edgecolor' : 20})
        plt.grid(False)
        plt.yticks([])
        plt.xticks([])
        plt.text(0.02, 0.5, 'total: %s / mean: %.2f / std: %.2f / min: %s / max: %s'\
            % (len(sequence_lengths),
               numpy.mean(sequence_lengths), numpy.std(sequence_lengths),\
               min(sequence_lengths),\
               max(sequence_lengths)),\
            va = 'center', alpha = 0.8, size = 'x-large')
    
        if dest == None:
            dest = self.fasta_file_path

        try:
            plt.savefig(dest + '.tiff')
        except:
            plt.savefig(dest + '.png')
    
        try:
            plt.show()
        except:
            pass
    
        return
    

if __name__ == '__main__':
    fasta = SequenceSource(sys.argv[1])
    output = FastaOutput(sys.argv[1] + 'out.fa')
    #fasta.visualize_sequence_length_distribution(title = sys.argv[2] if len(sys.argv) == 3 else 'None')
    while fasta.next():
        output.store(fasta)
    output.close()
