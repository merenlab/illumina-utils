# -*- coding: utf-8 -*-
#
# Copyright (C) 2013, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

import os
import sys
import gzip
import stat
import math
import numpy
import pickle
import textwrap

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.cm as cm
except:
    print("No matplotlib found")


import IlluminaUtils.lib.fastalib as u


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = remove_spaces(e)
        return
    def __str__(self):
        error_type = 'Config Error'

        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        error_lines = ['\033[0;30m\033[46m%s%s\033[0m' % (l, ' ' * (max_len - len(l)))\
                                         for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s: %s' % (error_type, error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(error_type) + 2), error_line))

        return '\n'.join(error_message)


def remove_spaces(text):
    while 1:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


conv_dict = {'A': 'T',
             'T': 'A',
             'C': 'G',
             'G': 'C',
             'N': 'N',
             '.': 'N'}


def is_file_exists(file_path):
    if not file_path:
        raise ConfigError("No input file is declared...")
    if not os.path.exists(file_path):
        raise ConfigError("No such file: '%s' :/" % file_path)
    return True


def is_output_file_writable(file_path):
    if not file_path:
        raise ConfigError("No output file is declared...")
    if not os.access(os.path.dirname(os.path.abspath(file_path)), os.W_OK):
        raise ConfigError("You do not have permission to generate the output file '%s'" % file_path)
    return True


def is_file_tab_delimited(file_path):
    is_file_exists(file_path)
    f = open(file_path)
    line = f.readline()
    if len(line.split('\t')) == 1:
        raise ConfigError("File '%s' does not seem to have TAB characters.\
                            Did you export this file on MAC using EXCEL? :(" % file_path)

    f.seek(0)
    if len(set([len(line.split('\t')) for line in f.readlines()])) != 1:
        raise ConfigError("Not all lines in the file '%s' have equal number of fields..." % file_path)

    f.close()
    return True


def get_TAB_delimited_file_as_dictionary(file_path, expected_fields = None, dict_to_append = None, column_names = None,\
                                        column_mapping = None, indexing_field = 0, assign_none_for_missing = False):
    is_file_exists(file_path)
    is_file_tab_delimited(file_path)

    f = open(file_path, 'rU')

    if column_names:
        columns = column_names
        num_fields = len(f.readline().strip('\n').split('\t'))

        if num_fields != len(columns):
            raise  ConfigError("Number of column names declared (%d) differs from the number of columns\
                                 found (%d) in the matrix ('%s') :/" % (len(columns), num_fields, file_path))
        f.seek(0)
    else:
        columns = [c.strip() for c in f.readline().strip('\n').split('\t')]

    if expected_fields:
        for field in expected_fields:
            if field not in columns:
                raise ConfigError("The file '%s' does not contain the right type of header. It was expected\
                                    to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(columns[1:])))

    d = {}

    for line in f.readlines():
        line_fields = line.strip('\n').strip('\r').split('\t')

        if column_mapping:
            updated_line_fields = []
            for i in range(0, len(line_fields)):
                try:
                    updated_line_fields.append(column_mapping[i](line_fields[i]))
                except NameError:
                    raise ConfigError("Mapping function '%s' did not work on value '%s'. These functions can be native\
                                        Python functions, such as 'str', 'int', or 'float', or anonymous functions\
                                        defined using lambda notation." % (column_mapping[i], line_fields[i]))
                except TypeError:
                    raise ConfigError("Mapping function '%s' does not seem to be a proper Python function :/" % column_mapping[i])
                except ValueError:
                    raise ConfigError("Mapping funciton '%s' did not like the value '%s' in column number %d\
                                        of the matrix :/" % (column_mapping[i], line_fields[i], i + 1))
            line_fields = updated_line_fields

        entry_name = line_fields[indexing_field]

        d[entry_name] = {}

        e = d[line_fields[indexing_field]]
        for i in range(0, len(columns)):
            if i == indexing_field:
                continue
            e[columns[i]] = line_fields[i]

    # we have the dict, but we will not return it the way it is if its supposed to be appended to an
    # already existing dictionary.
    if dict_to_append:
        # we don't want to through keys in d each time we want to add stuff to 'dict_to_append', so we keep keys we
        # find in the first item in the dict in another variable. this is potentially very dangerous if not every
        # item in 'd' has identical set of keys.
        keys = list(d.values())[0].keys()

        for entry in dict_to_append:
            if entry not in d:
                # so dict to append is missing a key that is in the dict to be appended. if the user did not
                # ask us to add None for these entries via none_for_missing, we are going to make a noise,
                # otherwise we will tolerate it.
                if not assign_none_for_missing:
                    raise ConfigError("Appending entries to the already existing dictionary from file '%s' failed\
                                        as the entry %s does not appear to be in the file." % (file_path, entry))
                else:
                    for key in keys:
                        dict_to_append[entry][key] = None
            else:
                for key in keys:
                    dict_to_append[entry][key] = d[entry][key]

        return dict_to_append

    return d


def colorize(sequence):
    Green = lambda s: '\033[30m\033[42m' + s + '' + '\033[0m'
    return Green(sequence)


def store_cPickle_obj(obj, output_file_path):
    f = gzip.open(output_file_path, 'wb')
    pickle.dump(obj, f)
    f.close()
    return True


def load_cPickle_obj(obj_file_path):
    f = gzip.open(obj_file_path, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj


def reverse_complement(seq):
    return ''.join(reversed([conv_dict[n] for n in seq]))


def reverse(seq):
    return ''.join(reversed(seq))


def complement(seq):
    return ''.join([conv_dict[n] for n in seq])



class ReadIDTracker:
    def __init__(self):
        self.ids = {}
        self.fates = set([])

    def update(self, pair_1, pair_2 = None, fate = 'unknown'):
        if fate not in self.fates:
            self.fates.add(fate)
            self.ids[fate] = set([])

        if pair_2:
            self.ids[fate].add((pair_1.entry.header_line, pair_2.entry.header_line),)
        else:
            self.ids[fate].add((pair_1.entry.header_line),)

    def store(self, output_file_path):
        store_cPickle_obj(self.ids, output_file_path)

class QualityScoresHandler:
    def __init__(self):
        self.data = {}
        self.entry_types = []
        self.finalized = False

    def update(self, pair_1, pair_2 = None, entry_type = 'default'):
        if entry_type not in self.entry_types:
            self.entry_types.append(entry_type)
            self.data[entry_type] = {}

        tile_number = pair_1.entry.tile_number
        tiles_dict = self.data[entry_type]

        q1 = pair_1.entry.process_Q_list()
        if pair_2:
            q2 = pair_2.entry.process_Q_list()

        if '1' not in tiles_dict:
            tiles_dict['1'] = {}
        if '2' not in tiles_dict:
            tiles_dict['2'] = {}

        if tile_number not in tiles_dict['1']:
            tiles_dict['1'][tile_number] = {'mean': [0] * len(q1), 'std': [0] * len(q1), 'count': [0] * len(q1)}
        if tile_number not in tiles_dict['2']:
            tiles_dict['2'][tile_number] = {'mean': [0] * len(q2), 'std': [0] * len(q2), 'count': [0] * len(q2)}

        tile_for_p1 = tiles_dict['1'][tile_number]
        if pair_2:
            tile_for_p2 = tiles_dict['2'][tile_number]

        # take care of the length variation:
        if len(q1) > len(tile_for_p1['mean']):
            diff = len(q1) - len(tile_for_p1['mean'])
            for attr in ['mean', 'std', 'count']:
                for d in range(0, diff):
                    tile_for_p1[attr].append(0)
        if pair_2 and len(q2) > len(tile_for_p2['mean']):
            diff = len(q2) - len(tile_for_p2['mean'])
            for attr in ['mean', 'std', 'count']:
                for d in range(0, diff):
                    tile_for_p2[attr].append(0)

        # update the tile
        for i in range(0, len(q1)):
            tile_for_p1['mean'][i] += q1[i]
            tile_for_p1['count'][i] += 1
        if pair_2:
            for i in range(0, len(q2)):
                tile_for_p2['mean'][i] += q2[i]
                tile_for_p2['count'][i] += 1

    def finalize(self):
        for entry_type in self.entry_types:
            for pair in ['1', '2']:
                for tile_id in self.data[entry_type][pair]:
                    tile = self.data[entry_type][pair][tile_id]
                    for i in range(0, len(tile['mean'])):
                        tile['mean'][i] = tile['mean'][i] * 1.0 / tile['count'][i]
        self.finalized = True

    def store_dict(self, output_file_path):
        if not self.finalized:
            self.finalize()

        store_cPickle_obj(self.data, output_file_path)


def quick_write(fp, header, sequence, qual):
    fp.write('@%s\n' % header)
    fp.write('%s\n' % sequence)
    fp.write('+\n')
    fp.write('%s\n' % qual)


def compute_plot_dict_from_tiles_dict(tiles_dict, plot_dict = {'1': {}, '2': {}}):
    sequence_length = max([len(t) for t in list(tiles_dict.values())[0].values()])

    for pair_no in ['1', '2']:
        for tile_no in tiles_dict[pair_no]:
            if tile_no not in plot_dict[pair_no]:
                plot_dict[pair_no][tile_no] = {'mean': [], 'count': [], 'std': []}

    for pair_no in ['1', '2']:
        for tile_no in tiles_dict[pair_no]:
            for i in range(0, sequence_length):
                plot_dict[pair_no][tile_no]['mean'].append(numpy.mean(tiles_dict[pair_no][tile_no][i]))
                plot_dict[pair_no][tile_no]['std'].append(numpy.std(tiles_dict[pair_no][tile_no][i]))
                plot_dict[pair_no][tile_no]['count'].append(len(tiles_dict[pair_no][tile_no][i]))

    return plot_dict

def visualize_sequence_length_distribution(fasta_file_path, dest, title, max_seq_len = None, xtickstep = None, ytickstep = None):
    sequence_lengths = []

    fasta = u.SequenceSource(fasta_file_path)

    while next(fasta):
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
    plt.fill_between(list(range(0, max_seq_len + 1)), seq_len_distribution, y2 = 0, color = 'black', alpha = 0.15)
    plt.ylabel('number of sequences')
    plt.xlabel('sequence length')

    if xtickstep == None:
        xtickstep = (max_seq_len / 50) or 1

    if ytickstep == None:
        ytickstep = max(seq_len_distribution) / 20 or 1

    plt.xticks(list(range(xtickstep, max_seq_len + 1, xtickstep)), rotation=90, size='xx-small')
    plt.yticks(list(range(0, max(seq_len_distribution) + 1, ytickstep)), size='xx-small')
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
        plt.savefig(dest + '.pdf')
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


def visualize_qual_stats_dict(D, dest, title, split_tiles = False, num_columns_to_show = 12, min_rows_to_show = 6):
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
    there are mean, std and count entries, in each of those, there are X items,
    where nth item corresponds to the mean of all values at nth location of the
    given tile.

    dest is destination file to save the output.

    title is the title to put on the figure.

    """

    # lets find out how many cycles were there. it is going to be about 101 for
    # hiseq runs, and 251 in miseq runs, but these values may change from run to
    # run. although all lanes are expected to have the identical number of cycles
    # the following code makes sure that the number_of_cycles variable holds the
    # longest one if there is a variation between the number of cycles between
    # lanes
    number_of_cycles = 0
    for pair in ['1', '2']:
        if pair not in D:
            continue
        for tile in D[pair]:
            if len(D[pair][tile]['mean']) > number_of_cycles:
                number_of_cycles = len(D[pair][tile]['mean'])


    def get_max_count(D, p = '1', tile = None):
        if tile == None:
            return float(max([D[p][x]['count'][0] for x in list(D[p].keys())]))
        else:
            return float(max(D[p][tile]['count']))


    def get_num_tiles(D, p = '1'):
        return len(list(D[p].keys()))

    num_tiles = get_num_tiles(D)

    if split_tiles:
        num_rows_to_show = int((math.ceil(num_tiles / num_columns_to_show)) or 1) * 2
    else:
        num_rows_to_show = int((math.ceil(num_tiles / num_columns_to_show)) or 1)

    if num_rows_to_show < min_rows_to_show:
        num_rows_to_show = min_rows_to_show

    fig = plt.figure(figsize = (3 * num_columns_to_show, 2 * num_rows_to_show))
    gs = Gs(num_rows_to_show, num_columns_to_show)

    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)


    subplots = {}
    colors = cm.get_cmap('RdYlGn', lut=256)

    tiles = sorted(D['1'].keys())

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

                plt.xticks(list(range(number_of_cycles / 10, number_of_cycles, number_of_cycles / 10)), rotation=90, size='xx-small')
                plt.ylim(ymin = 0, ymax = 42)
                plt.xlim(xmin = 0, xmax = number_of_cycles - 1)
                if _pair == '1':
                    plt.yticks(list(range(5, 41, 5)), size='xx-small')
                else:
                    plt.yticks(list(range(5, 41, 5)), [])

                if tile in D[_pair]:
                    plt.fill_between(list(range(0, number_of_cycles)), [42 for _ in range(0, number_of_cycles)], y2 = 0, color = colors(D[_pair][tile]['count'][0] / m[_pair]), alpha = 0.2)
                    subplots[tile][_pair].plot(D[_pair][tile]['mean'], color = _color, lw = 2)

                    read_number_percent_dropdown = [42 * (x / get_max_count(D, _pair, tile)) for x in D[_pair][tile]['count']]
                    if not len(set(read_number_percent_dropdown)) <= 1:
                        plt.fill_between(list(range(0, number_of_cycles)), read_number_percent_dropdown, y2 = 0, color = 'black', alpha = 0.08)

                    plt.text(5, 2.5, '%s/%s :: %s' % (tile, _pair, big_number_pretty_print(int(get_max_count(D, _pair, tile)))), alpha=0.8, size = 'x-small')
                else:
                    plt.text(5, 2.5, '%s/%s :: 0' % (tile, _pair), alpha=0.8, size = 'x-small')
                    plt.fill_between(list(range(0, number_of_cycles)), [42 for _ in range(0, number_of_cycles)], y2 = 0, color = colors(0 / m[_pair]), alpha = 0.2)

        Total = lambda p: big_number_pretty_print(int(sum([D[p][x]['count'][0] for x in D[p]])))

        plt.figtext(0.5, 0.97, '%s (p1: %s; p2: %s)' % (title, Total('1'), Total('2')), weight = 'black', size = 'xx-large', ha = 'center')

    else:
        m = get_max_count(D)

        for i in range(0, len(tiles)):
            tile = tiles[i]

            subplots[tile] = plt.subplot(next(gs))
            plt.grid(True)

            plt.subplots_adjust(left=0.02, bottom = 0.03, top = 0.95, right = 0.98)

            plt.xticks(list(range(number_of_cycles / 10, number_of_cycles, number_of_cycles / 10)), rotation=90, size='xx-small')
            plt.ylim(ymin = 0, ymax = 42)
            plt.xlim(xmin = 0, xmax = number_of_cycles - 1)
            plt.yticks(list(range(5, 41, 5)), size='xx-small')

            if tile in D['1']:
                plt.fill_between(list(range(0, number_of_cycles)), [42 for _ in range(0, number_of_cycles)], y2 = 0, color = colors(D['1'][tile]['count'][0] / m), alpha = 0.2)
                subplots[tile].plot(D['1'][tile]['mean'], color = 'orange', lw = 2)

                read_number_percent_dropdown = [42 * (x / get_max_count(D, tile = tile)) for x in D['1'][tile]['count']]
                if not len(set(read_number_percent_dropdown)) <= 1:
                    if len(read_number_percent_dropdown) < number_of_cycles:
                        read_number_percent_dropdown += [0] * (number_of_cycles - len(read_number_percent_dropdown))
                    plt.fill_between(list(range(0, number_of_cycles)), read_number_percent_dropdown, y2 = 0, color = 'black', alpha = 0.08)

                plt.text(5, 2.5, '%s :: %s' % (tile, big_number_pretty_print(int(get_max_count(D, tile = tile)))), alpha=0.5)
            else:
                plt.text(5, 2.5, '%s :: 0' % tile, alpha=0.5)
                plt.fill_between(list(range(0, number_of_cycles)), [42 for _ in range(0, number_of_cycles)], y2 = 0, color = colors(0 / m), alpha = 0.2)
            if tile in D['2']:
                subplots[tile].plot(D['2'][tile]['mean'], color = 'purple', lw = 2)

        plt.figtext(0.5, 0.97, '%s (%s)' % (title, big_number_pretty_print(int(sum([D['1'][x]['count'][0] for x in D['1']])))), weight = 'black', size = 'xx-large', ha = 'center')

    plt.savefig(dest + '.png')


def visualize_qual_stats_dict_single(D, dest, title):
    """
    same as visualize_qual_stats_dict, but puts all tiles together.
    """

    # first find out how many cycles were there. it is going to be about 101 for
    # hiseq runs, and 251 in miseq runs, but these values may change from run to
    # run. although all lanes are expected to have the identical number of cycles
    # the following code makes sure that the number_of_cycles variable holds the
    # longest one if there is a variation between the number of cycles between
    # lanes
    number_of_cycles = 0
    for pair in ['1', '2']:
        if pair not in D:
            continue
        for tile in D[pair]:
            if len(D[pair][tile]['mean']) > number_of_cycles:
                number_of_cycles = len(D[pair][tile]['mean'])


    fig = plt.figure(figsize = (12, 8))

    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)

    all_tiles = {'1': {'mean': [0] * number_of_cycles, 'count': [0] * number_of_cycles},
                 '2': {'mean': [0] * number_of_cycles, 'count': [0] * number_of_cycles}
                }

    for i in range(0, number_of_cycles):
        means_p1 = []
        counts_p1 = []
        means_p2 = []
        counts_p2 = []
        for tile_id in D['1']:
            tile = D['1'][tile_id]
            means_p1.append(tile['mean'][i])
            counts_p1.append(tile['count'][i])
            if '2' in D and D['2']:
                tile = D['2'][tile_id]
                means_p2.append(tile['mean'][i])
                counts_p2.append(tile['count'][i])

        all_tiles['1']['mean'][i] = numpy.mean(means_p1)
        all_tiles['1']['count'][i] = sum(counts_p1)

        if '2' in D and D['2']:
            all_tiles['2']['mean'][i] = numpy.mean(means_p2)
            all_tiles['2']['count'][i] = sum(counts_p2)

    colors = cm.get_cmap('RdYlGn', lut=256)

    plt.grid(True)

    plt.subplots_adjust(left=0.02, bottom = 0.03, top = 0.95, right = 0.98)

    plt.xticks(list(range(number_of_cycles / 10, number_of_cycles, number_of_cycles / 10)), rotation=90, size='xx-small')
    plt.ylim(ymin = 0, ymax = 42)
    plt.xlim(xmin = 0, xmax = number_of_cycles - 1)
    plt.yticks(list(range(5, 41, 5)), size='xx-small')

    plt.fill_between(list(range(0, number_of_cycles)), [42 for _ in range(0, number_of_cycles)], y2 = 0, color = colors(0), alpha = 0.2)
    plt.plot(all_tiles['1']['mean'], color = 'orange', lw = 6)

    read_number_percent_dropdown = [42 * (x / all_tiles['1']['count'][0]) for x in all_tiles['1']['count']]
    if not len(set(read_number_percent_dropdown)) <= 1:
        plt.fill_between(list(range(0, number_of_cycles)), read_number_percent_dropdown, y2 = 0, color = 'black', alpha = 0.08)

        plt.text(5, 2.5, '%s' % big_number_pretty_print(all_tiles['1']['count'][0]), alpha=0.5)
    else:
        plt.text(5, 2.5, '%s' % big_number_pretty_print(all_tiles['1']['count'][0]), alpha=0.5)

    if '2' in all_tiles and all_tiles['2']:
        plt.plot(all_tiles['2']['mean'], color = 'purple', lw = 6)

    plt.figtext(0.5, 0.97, '%s' % (title), weight = 'black', size = 'xx-large', ha = 'center')

    try:
        plt.savefig(dest + '.tiff')
    except:
        plt.savefig(dest + '.png')

    return (all_tiles['1']['mean'], all_tiles['2']['mean'])


def populate_tiles_qual_dict_from_input(input_1, input_2, tiles_dict = {'1': {}, '2': {}}):
    while next(input_1):
        if input_1.p_available:
            input_1.print_percentage()

        q1 = input_1.entry.process_Q_list()

        if input_1.entry.tile_number not in tiles_dict['1']:
            tiles_dict['1'][input_1.entry.tile_number] = []
            for i in range(0, len(q1)):
                tiles_dict['1'][input_1.entry.tile_number].append([])

        for i in range(0, len(q1)):
            try:
                tiles_dict['1'][input_1.entry.tile_number][i].append(q1[i])
            except:
                tiles_dict['1'][input_1.entry.tile_number].append([])
                tiles_dict['1'][input_1.entry.tile_number][i].append(q1[i])


    sys.stderr.write('\n')

    while next(input_2):
        if input_2.p_available:
            input_2.print_percentage()

        if input_2.entry.tile_number not in tiles_dict['2']:
            tiles_dict['2'][input_2.entry.tile_number] = []
            for i in range(0, 101):
                tiles_dict['2'][input_2.entry.tile_number].append([])

        q2 = input_2.entry.process_Q_list()

        for i in range(0, len(q2)):
            try:
                tiles_dict['2'][input_2.entry.tile_number][i].append(q2[i])
            except:
                tiles_dict['2'][input_2.entry.tile_number].append([])

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
