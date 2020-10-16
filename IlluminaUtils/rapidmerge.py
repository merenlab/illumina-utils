#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2020, Samuel E. Miller
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import io
import os
import re
import sys
import gzip
import stat
import struct
import functools
import multiprocessing

from datetime import datetime
from multiprocessing import Process, Queue, Value

from IlluminaUtils.utils.terminal import Progress
from IlluminaUtils.lib.fastqlib import FastQSource, FastQEntry
from IlluminaUtils.utils.helperfunctions import (
    big_number_pretty_print, combine_files, conv_dict, is_file_exists, reverse_complement)

try:
    import Levenshtein
except:
    print('''
    ERROR: This script requires Levenshtein module, which seems to not installed on this machine.
    Here is a fast implementation of Levenshtein distance for Python:
        http://code.google.com/p/pylevenshtein/
    ''')
    sys.exit(0)


class Multiprocessor:

    def __init__(self):
        self.processes = []
        self.queue = Queue()
        return


    @staticmethod
    def _wrapper(func, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put(ret)
        return


    def run(self, func, *args, **kwargs):
        p = Process(target=self._wrapper, args=[func, self.queue, args, kwargs])
        self.processes.append(p)
        p.start()
        return


    def wait(self):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        for p in self.processes:
            p.join()
        return rets


class FASTQMerger:
    def __init__(
        self,
        input1_path,
        input2_path,
        input1_is_gzipped=False,
        input2_is_gzipped=False,
        ignore_deflines=False,
        output_dir='',
        output_file_name='output',
        r1_prefix_pattern='',
        r2_prefix_pattern='',
        report_r1_prefix=False,
        report_r2_prefix=False,
        max_p=0.3,
        max_num_mismatches=-1,
        min_overlap_size=16,
        min_qual_score=15,
        partial_overlap_only=True,
        retain_overlap_only=False,
        skip_suffix_trimming=False,
        ignore_Ns=False,
        enforce_Q30_check=False,
        dataset_index=None,
        total_dataset_count=None,
        project_name='Unknown_project',
        num_cores=1):

        is_file_exists(input1_path)
        is_file_exists(input2_path)
        self.input1_path = input1_path
        self.input2_path = input2_path
        self.input1_is_gzipped = input1_is_gzipped
        self.input2_is_gzipped = input2_is_gzipped
        self.ignore_deflines = ignore_deflines

        # Check the validity of the FASTQ files and get the lengths of trimmed reads.
        fastq_source = FastQSource(input1_path, compressed=input1_is_gzipped)
        fastq_source.next()
        first_seq1 = fastq_source.entry.sequence
        fastq_source.close()
        fastq_source = FastQSource(input2_path, compressed=input2_is_gzipped)
        fastq_source.next()
        first_seq2 = fastq_source.entry.sequence
        fastq_source.close()

        self.read1_length = len(first_seq1)
        self.read2_length = len(first_seq2)

        output_path_maker = lambda p: os.path.join(output_dir, output_file_name + p)
        self.merged_path = output_path_maker('_MERGED')
        self.merge_failed_path = output_path_maker('_FAILED')
        self.merge_failed_with_Ns_path = output_path_maker('_FAILED_WITH_Ns')

        self.r1_prefix_pattern = r1_prefix_pattern
        self.r2_prefix_pattern = r2_prefix_pattern
        self.r1_prefix_compiled = re.compile(r1_prefix_pattern)
        self.r2_prefix_compiled = re.compile(r2_prefix_pattern)

        if report_r1_prefix and r1_prefix_pattern == '':
            raise UserWarning(
                "In the absence of a read 1 prefix pattern, "
                "no file of read 1 prefix sequences from merged reads will be reported, "
                "despite a `--report_r1_prefix` flag.")
            report_r1_prefix = False
        if report_r2_prefix and r2_prefix_pattern == '':
            raise UserWarning(
                "In the absence of a read 2 prefix pattern, "
                "no file of read 2 prefix sequences from merged reads will be reported, "
                "despite a `--report_r2_prefix` flag.")
            report_r2_prefix = False
        self.report_r1_prefix = report_r1_prefix
        self.report_r2_prefix = report_r2_prefix
        self.r1_prefix_path = output_path_maker('_MERGED_R1_PREFIX') if report_r1_prefix else ''
        self.r2_prefix_path = output_path_maker('_MERGED_R2_PREFIX') if report_r2_prefix else ''

        self.max_p = max_p
        self.max_num_mismatches = max_num_mismatches
        self.min_overlap_size = min_overlap_size
        self.min_qual_score = min_qual_score
        self.partial_overlap_only = partial_overlap_only
        self.retain_overlap_only = retain_overlap_only
        self.skip_suffix_trimming = skip_suffix_trimming
        self.ignore_Ns = ignore_Ns

        self.enforce_Q30_check = enforce_Q30_check
        self.merge_failed_Q30_path = output_path_maker('_FAILED_Q30') if enforce_Q30_check else ''

        self.dataset_index = dataset_index
        self.total_dataset_count = total_dataset_count
        self.project_name = project_name

        if not 0 < num_cores <= multiprocessing.cpu_count():
            raise RuntimeError("\"%d\" is not a valid number of cores. "
                               "The number of cores must be between 1 and %d."
                               % (num_cores, multiprocessing.cpu_count()))
        self.num_cores = num_cores

        self.stats = FASTQMergerStats()
        self.stats_path = output_path_maker('_STATS')
        self.mismatches_breakdown_path = output_path_maker('_MISMATCHES_BREAKDOWN')

        return


    def run(self, merge_method='hamming'):

        # values for tracking progress
        num_lines_total = 0
        input1_file = (gzip.open(self.input1_path) if self.input1_is_gzipped
                       else open(self.input1_path))
        for line in input1_file:
            num_lines_total += 1
        input1_file.close()
        num_items_total = num_lines_total // 4
        num_items_processed = Value('i', 0)
        next_percentage = Value('i', 1)
        num_pairs_prefix_passed = Value('i', 0)
        num_merged_pairs = Value('i', 0)
        num_merged_pairs_zero_mismatches = Value('i', 0)

        # singlethreaded
        if self.num_cores == 1:
            self.stats.update(merge_reads_in_files(
                self.input1_path,
                self.input2_path,
                self.read1_length,
                self.read2_length,
                merge_method=merge_method,
                input1_is_gzipped=self.input1_is_gzipped,
                input2_is_gzipped=self.input2_is_gzipped,
                ignore_deflines=self.ignore_deflines,
                merged_path=self.merged_path,
                merge_failed_path=self.merge_failed_path,
                merge_failed_with_Ns_path=self.merge_failed_with_Ns_path,
                merge_failed_Q30_path=self.merge_failed_Q30_path,
                r1_prefix_compiled=self.r1_prefix_compiled,
                r2_prefix_compiled=self.r2_prefix_compiled,
                r1_prefix_path=self.r1_prefix_path,
                r2_prefix_path=self.r2_prefix_path,
                max_p=self.max_p,
                max_num_mismatches=self.max_num_mismatches,
                min_overlap_size=self.min_overlap_size,
                min_qual_score=self.min_qual_score,
                partial_overlap_only=self.partial_overlap_only,
                retain_overlap_only=self.retain_overlap_only,
                skip_suffix_trimming=self.skip_suffix_trimming,
                ignore_Ns=self.ignore_Ns,
                enforce_Q30_check=self.enforce_Q30_check,
                dataset_index=self.dataset_index,
                total_dataset_count=self.total_dataset_count,
                project_name=self.project_name,
                verbose=True,
                num_items_total=num_items_total,
                num_items_processed=num_items_processed,
                next_percentage=next_percentage,
                num_pairs_prefix_passed=num_pairs_prefix_passed,
                num_merged_pairs=num_merged_pairs,
                num_merged_pairs_zero_mismatches=num_merged_pairs_zero_mismatches))
            self.stats.write_stats(
                self.stats_path,
                self.max_p,
                self.max_num_mismatches,
                self.min_overlap_size,
                self.min_qual_score,
                self.ignore_Ns,
                self.enforce_Q30_check,
                self.partial_overlap_only,
                self.retain_overlap_only)
            self.stats.write_mismatches_breakdown_table(self.mismatches_breakdown_path)
            print()
            return

        # Prepare multiprocessing.
        # Find positions at which to chunk the input files.
        r1_start_positions, r1_end_positions, r1_end_strings, r2_start_positions = self.find_fastq_chunk_starts()
        # Create temporary output files.
        time_str = datetime.now().isoformat(timespec='seconds').replace('-', '').replace(':', '')
        temp_merged_paths = [
            ("%s_TEMP_%d_%s"
             % (self.merged_path, chunk_index, time_str))
            for chunk_index in range(1, self.num_cores + 1)]
        temp_merge_failed_paths = [
            ("%s_TEMP_%d_%s"
             % (self.merge_failed_path, chunk_index, time_str))
            for chunk_index in range(1, self.num_cores + 1)]
        temp_merge_failed_with_Ns_paths = [
            ("%s_TEMP_%d_%s"
             % (self.merge_failed_with_Ns_path, chunk_index, time_str))
            for chunk_index in range(1, self.num_cores + 1)]
        temp_merge_failed_Q30_paths = [
            ("%s_TEMP_%d_%s"
             % (self.merge_failed_Q30_path, chunk_index, time_str))
            if self.merge_failed_Q30_path else ''
            for chunk_index in range(1, self.num_cores + 1)]
        temp_r1_prefix_paths = [
            ("%s_TEMP_%d_%s"
             % (self.r1_prefix_path, chunk_index, time_str))
            if self.r1_prefix_path else ''
            for chunk_index in range(1, self.num_cores + 1)]
        temp_r2_prefix_paths = [
            ("%s_TEMP_%d_%s"
             % (self.r2_prefix_path, chunk_index, time_str))
            if self.r2_prefix_path else ''
            for chunk_index in range(1, self.num_cores + 1)]

        mp = Multiprocessor()
        for (temp_merged_path,
             temp_merge_failed_path,
             temp_merge_failed_with_Ns_path,
             temp_merge_failed_Q30_path,
             temp_r1_prefix_path,
             temp_r2_prefix_path,
             r1_start_position,
             r1_end_position,
             r1_end_string,
             r2_start_position) in zip(temp_merged_paths,
                                       temp_merge_failed_paths,
                                       temp_merge_failed_with_Ns_paths,
                                       temp_merge_failed_Q30_paths,
                                       temp_r1_prefix_paths,
                                       temp_r2_prefix_paths,
                                       r1_start_positions,
                                       r1_end_positions,
                                       r1_end_strings,
                                       r2_start_positions):
            mp.run(
                merge_reads_in_files,
                *(
                    self.input1_path,
                    self.input2_path,
                    self.read1_length,
                    self.read2_length),
                **{
                    'merge_method': merge_method,
                    'input1_is_gzipped': self.input1_is_gzipped,
                    'input2_is_gzipped': self.input2_is_gzipped,
                    'ignore_deflines': self.ignore_deflines,
                    'merged_path': temp_merged_path,
                    'merge_failed_path': temp_merge_failed_path,
                    'merge_failed_with_Ns_path': temp_merge_failed_with_Ns_path,
                    'merge_failed_Q30_path': temp_merge_failed_Q30_path,
                    'r1_prefix_compiled': self.r1_prefix_compiled,
                    'r2_prefix_compiled': self.r2_prefix_compiled,
                    'r1_prefix_path': temp_r1_prefix_path,
                    'r2_prefix_path': temp_r2_prefix_path,
                    'max_p': self.max_p,
                    'max_num_mismatches': self.max_num_mismatches,
                    'min_overlap_size': self.min_overlap_size,
                    'min_qual_score': self.min_qual_score,
                    'partial_overlap_only': self.partial_overlap_only,
                    'retain_overlap_only': self.retain_overlap_only,
                    'skip_suffix_trimming': self.skip_suffix_trimming,
                    'ignore_Ns': self.ignore_Ns,
                    'enforce_Q30_check': self.enforce_Q30_check,
                    'dataset_index': self.dataset_index,
                    'total_dataset_count': self.total_dataset_count,
                    'project_name': self.project_name,
                    'verbose': True,
                    'num_items_total': num_items_total,
                    'num_items_processed': num_items_processed,
                    'next_percentage': next_percentage,
                    'num_pairs_prefix_passed': num_pairs_prefix_passed,
                    'num_merged_pairs': num_merged_pairs,
                    'num_merged_pairs_zero_mismatches': num_merged_pairs_zero_mismatches,
                    'r1_start_position': r1_start_position,
                    'r1_end_position': r1_end_position,
                    'r1_end_string': r1_end_string,
                    'r2_start_position': r2_start_position})
        chunk_stats = mp.wait()
        print()

        # Delete temp files after combining them.
        combine_files(temp_merged_paths, self.merged_path)
        for path in temp_merged_paths:
            if os.path.exists(path):
                os.remove(path)
        combine_files(temp_merge_failed_paths, self.merge_failed_path)
        for path in temp_merge_failed_paths:
            if os.path.exists(path):
                os.remove(path)
        combine_files(temp_merge_failed_with_Ns_paths, self.merge_failed_with_Ns_path)
        for path in temp_merge_failed_with_Ns_paths:
            if os.path.exists(path):
                os.remove(path)
        if self.merge_failed_Q30_path:
            combine_files(temp_merge_failed_Q30_paths, self.merge_failed_Q30_path)
            for path in temp_merge_failed_Q30_paths:
                if os.path.exists(path):
                    os.remove(path)
        if self.r1_prefix_path:
            combine_files(temp_r1_prefix_paths, self.r1_prefix_path)
            for path in temp_r1_prefix_paths:
                if os.path.exists(path):
                    os.remove(path)
        if self.r2_prefix_path:
            combine_files(temp_r2_prefix_paths, self.r2_prefix_path)
            for path in temp_r2_prefix_paths:
                if os.path.exists(path):
                    os.remove(path)

        for s in chunk_stats:
            if s:
                self.stats.update(s)
        self.stats.write_stats(
            self.stats_path,
            self.max_p,
            self.max_num_mismatches,
            self.min_overlap_size,
            self.min_qual_score,
            self.ignore_Ns,
            self.enforce_Q30_check,
            self.partial_overlap_only,
            self.retain_overlap_only)
        self.stats.write_mismatches_breakdown_table(self.mismatches_breakdown_path)
        return


    def find_fastq_chunk_starts(self):

        if self.input1_is_gzipped:
            with gzip.open(self.input1_path, 'rt') as f:
                uncompressed_file_size = f.seek(0, io.SEEK_END)
        else:
            uncompressed_file_size = os.stat(self.input1_path)[stat.ST_SIZE]
        chunk_size = uncompressed_file_size // self.num_cores

        fastq_file1 = gzip.open(self.input1_path, 'rt') if self.input1_is_gzipped else open(self.input1_path)
        r1_start_positions = []
        r1_end_strings = []
        position = 0
        prev_position = 0
        for chunk in range(self.num_cores):
            fastq_file1.seek(position)
            while True:
                line_end = fastq_file1.readline()
                # Checking for empty string must come before checking for '@'
                if line_end == '':
                    # If EOF, append -1 rather than the last position.
                    r1_start_positions.append(-1)
                    break
                prev_position = position
                position = fastq_file1.tell()
                if line_end[0] == '@':
                    r1_start_positions.append(prev_position)
                    position += chunk_size
                    break
            if chunk > 0:
                r1_end_strings.append(line_end.rstrip())
        r1_end_strings.append('')
        r1_end_positions = [p for p in r1_start_positions[1:]] + [-1]

        if self.read1_length != self.read2_length:
            r2_start_positions = [0]
            fastq_file2 = gzip.open(self.input2_path, 'rt') if self.input2_is_gzipped else open(self.input2_path)
            r1_end_string_iter = iter(r1_end_strings)
            e = next(r1_end_string_iter)
            if ' 1:' in e:
                e = e.replace(' 1:', ' 2:') + '\n'
            else:
                e = e.replace(' 2:', ' 1:') + '\n'
            line = True
            while line:
                line = fastq_file2.readline()
                if line == e:
                    position = fastq_file2.tell() - len(line)
                    r2_start_positions.append(position)
                    e = next(r1_end_string_iter)
                    if ' 1:' in e:
                        e = e.replace(' 1:', ' 2:') + '\n'
                    else:
                        e = e.replace(' 2:', ' 1:') + '\n'
            fastq_file2.close()
        else:
            r2_start_positions = r1_start_positions

        fastq_file1.close()

        return r1_start_positions, r1_end_positions, r1_end_strings, r2_start_positions


class FASTQMergerStats:
    def __init__(self):
        self.actual_number_of_pairs = 0
        self.pairs_eliminated_due_to_Ns = 0
        self.pairs_eliminated_due_to_P = 0
        self.pairs_eliminated_due_to_max_num_mismatches = 0
        self.pairs_eliminated_due_to_Q30 = 0
        self.pairs_eliminated_due_to_Min_Overlap = 0
        self.prefix_failed_in_pair_1_total = 0
        self.prefix_failed_in_pair_2_total = 0
        self.prefix_failed_in_both_pairs_total = 0
        self.passed_prefix_total = 0
        self.failed_prefix_total = 0
        self.merge_failed_total = 0
        self.merge_passed_total = 0
        self.merging_done_on_complete_overlap = 0
        self.num_mismatches_breakdown = {}
        self.num_mismatches_breakdown['merge passed'] = {}
        self.num_mismatches_breakdown['merge passed'][0] = 0
        self.total_number_of_mismatches = 0
        self.num_mismatches_recovered_from_read_1 = 0
        self.num_mismatches_recovered_from_read_2 = 0
        self.num_mismatches_replaced_with_N = 0
        self.num_Q30_fails_in_read_1 = 0
        self.num_Q30_fails_in_read_2 = 0
        self.num_Q30_fails_in_both = 0
        return


    def update(self, stats):
        self.actual_number_of_pairs += stats.actual_number_of_pairs
        self.pairs_eliminated_due_to_Ns += stats.pairs_eliminated_due_to_Ns
        self.pairs_eliminated_due_to_P += stats.pairs_eliminated_due_to_P
        self.pairs_eliminated_due_to_max_num_mismatches += (
            stats.pairs_eliminated_due_to_max_num_mismatches)
        self.pairs_eliminated_due_to_Q30 += stats.pairs_eliminated_due_to_Q30
        self.pairs_eliminated_due_to_Min_Overlap += stats.pairs_eliminated_due_to_Min_Overlap
        self.prefix_failed_in_pair_1_total += stats.prefix_failed_in_pair_1_total
        self.prefix_failed_in_pair_2_total += stats.prefix_failed_in_pair_2_total
        self.prefix_failed_in_both_pairs_total += stats.prefix_failed_in_both_pairs_total
        self.passed_prefix_total += stats.passed_prefix_total
        self.failed_prefix_total += stats.failed_prefix_total
        self.merge_failed_total += stats.merge_failed_total
        self.merge_passed_total += stats.merge_passed_total
        self.merging_done_on_complete_overlap += stats.merging_done_on_complete_overlap
        merge_passed_breakdown = self.num_mismatches_breakdown['merge passed']
        for i, j in stats.num_mismatches_breakdown['merge passed'].items():
            if i in merge_passed_breakdown:
                merge_passed_breakdown[i] += j
            else:
                merge_passed_breakdown[i] = j
        self.total_number_of_mismatches += stats.total_number_of_mismatches
        self.num_mismatches_recovered_from_read_1 += stats.num_mismatches_recovered_from_read_1
        self.num_mismatches_recovered_from_read_2 += stats.num_mismatches_recovered_from_read_2
        self.num_mismatches_replaced_with_N += stats.num_mismatches_replaced_with_N
        self.num_Q30_fails_in_read_1 += stats.num_Q30_fails_in_read_1
        self.num_Q30_fails_in_read_2 += stats.num_Q30_fails_in_read_2
        self.num_Q30_fails_in_both += stats.num_Q30_fails_in_both
        return


    def s_line(self, label, value, padding=55):
        return '%s%s\t%s\n' % (label, ' ' + '.' * (padding - len(label)), value)


    def write_stats(
        self,
        output_file_path,
        max_p,
        max_num_mismatches,
        min_overlap_size,
        min_qual_score,
        ignore_Ns,
        enforce_Q30_check,
        partial_overlap_only,
        retain_overlap_only):

        stats = open(output_file_path, 'w')

        stats.write(self.s_line('Number of pairs analyzed', '%d' % self.actual_number_of_pairs))
        stats.write(self.s_line(
            'Prefix failed in read 1', '%d' % self.prefix_failed_in_pair_1_total))
        stats.write(
            self.s_line('Prefix failed in read 2', '%d' % self.prefix_failed_in_pair_2_total))
        stats.write(
            self.s_line('Prefix failed in both', '%d' % self.prefix_failed_in_both_pairs_total))
        stats.write(self.s_line('Passed prefix total', '%d' % self.passed_prefix_total))
        stats.write(self.s_line('Failed prefix total', '%d' % self.failed_prefix_total))
        stats.write(self.s_line('Merged total', '%d' % self.merge_passed_total))
        stats.write(self.s_line(
            'Complete overlap forced total', '%d' % self.merging_done_on_complete_overlap))
        stats.write(self.s_line('Merge failed total', '%d' % self.merge_failed_total))
        stats.write(self.s_line('Merge discarded due to P', '%d' % self.pairs_eliminated_due_to_P))
        stats.write(self.s_line(
            'Merge discarded due to max num mismatches',
            '%d' % self.pairs_eliminated_due_to_max_num_mismatches))
        stats.write(self.s_line(
            'Merge discarded due to Ns', '%d' % self.pairs_eliminated_due_to_Ns))
        stats.write(self.s_line('Merge discarded due to Q30', '%d' % self.pairs_eliminated_due_to_Q30))
        stats.write(self.s_line(
            'Pairs discarded due to min expected overlap',
            '%d' % self.pairs_eliminated_due_to_Min_Overlap))
        stats.write(self.s_line(
            'Num mismatches found in merged reads', '%d' % self.total_number_of_mismatches))
        stats.write(self.s_line(
            'Mismatches recovered from read 1', '%d' % self.num_mismatches_recovered_from_read_1))
        stats.write(self.s_line(
            'Mismatches recovered from read 2', '%d' % self.num_mismatches_recovered_from_read_2))
        stats.write(self.s_line(
            'Mismatches replaced with N', '%d' % self.num_mismatches_replaced_with_N))
        stats.write('\n\nMismatches breakdown in final merged reads:\n\n')

        for i in sorted(self.num_mismatches_breakdown['merge passed']):
            stats.write('%d\t%d\n' % (i, self.num_mismatches_breakdown['merge passed'][i]))

        stats.write('\n\n')
        stats.write(self.s_line('Command line', '%s' % ' '.join(sys.argv)))
        stats.write(self.s_line('Work directory', '%s' % os.getcwd()))
        stats.write(self.s_line('"p" value', '%f' % max_p))
        stats.write(self.s_line(
            'Max num mismatches at the overlapped region',
            '%s' % (str(max_num_mismatches) if max_num_mismatches >= 0 else 'None')))
        stats.write(self.s_line('Min overlap size', '%s' % min_overlap_size))
        stats.write(self.s_line('Min Q-score for mismatches', '%s' % min_qual_score))
        stats.write(self.s_line('Ns ignored?', '%s' % ignore_Ns))
        stats.write(self.s_line('Q30 enforced?', '%s' % enforce_Q30_check))
        stats.write(self.s_line(
            'Run with stringent flag on?', '%s' % (not partial_overlap_only)))
        stats.write(self.s_line(
            'Requested only the overlapping part to be retained?', '%s' % retain_overlap_only))

        stats.close()
        return


    def write_mismatches_breakdown_table(self, output_file_path):
        num_mismatches_breakdown_table = open(output_file_path, 'w')
        categories = list(self.num_mismatches_breakdown.keys())
        num_mismatches_breakdown_table.write('%s\t%s\n' % ('num_mismatch', '\t'.join(categories)))

        for i in range(0, max([max(self.num_mismatches_breakdown[x].keys())
                               for x in self.num_mismatches_breakdown])):
            mismatch_counts = []
            for category in categories:
                count = 0
                if i in self.num_mismatches_breakdown[category]:
                    count = self.num_mismatches_breakdown[category][i]

                mismatch_counts.append(count)

            num_mismatches_breakdown_table.write(
                '%d\t%s\n' % (i, '\t'.join([str(x) for x in mismatch_counts])))

        num_mismatches_breakdown_table.close()
        return


    def record_num_mismatches(self, number_of_mismatches, merging_result='merged'):
        if merging_result not in self.num_mismatches_breakdown:
            self.num_mismatches_breakdown[merging_result] = {}

        if number_of_mismatches not in self.num_mismatches_breakdown[merging_result]:
            self.num_mismatches_breakdown[merging_result][number_of_mismatches] = 1
        else:
            self.num_mismatches_breakdown[merging_result][number_of_mismatches] += 1
        return


    def process_recovery_dict(self, recovery_dict):
        self.total_number_of_mismatches += sum(recovery_dict.values())
        self.num_mismatches_recovered_from_read_1 += recovery_dict['r1']
        self.num_mismatches_recovered_from_read_2 += recovery_dict['r2']
        self.num_mismatches_replaced_with_N += recovery_dict['none']
        return


def merge_reads_in_files(
    input1_path,
    input2_path,
    read1_length,
    read2_length,
    merge_method='hamming',
    progress_info=None,
    input1_is_gzipped=False,
    input2_is_gzipped=False,
    ignore_deflines=False,
    merged_path='output_MERGED',
    merge_failed_path='output_FAILED',
    merge_failed_with_Ns_path='output_FAILED_WITH_Ns',
    merge_failed_Q30_path='',
    r1_prefix_compiled=None,
    r2_prefix_compiled=None,
    r1_prefix_path='',
    r2_prefix_path='',
    max_p=0.3,
    max_num_mismatches=-1,
    min_overlap_size=16,
    min_qual_score=15,
    partial_overlap_only=True,
    retain_overlap_only=False,
    skip_suffix_trimming=False,
    ignore_Ns=False,
    enforce_Q30_check=False,
    dataset_index=None,
    total_dataset_count=None,
    project_name='Unknown_project',
    verbose=True,
    num_items_total=None,
    num_items_processed=None,
    next_percentage=None,
    num_pairs_prefix_passed=None,
    num_merged_pairs=None,
    num_merged_pairs_zero_mismatches=None,
    r1_start_position=0,
    r1_end_position=-1,
    r1_end_string='',
    r2_start_position=0):

    if merge_method == 'hamming':
        partial_overlap_merge_function = functools.partial(merge_by_distance_metric, metric=Levenshtein.hamming)
    elif merge_method == 'levenshtein':
        partial_overlap_merge_function = functools.partial(merge_by_distance_metric, metric=Levenshtein.distance)
    elif merge_method == 'exact':
        partial_overlap_merge_function = merge_with_zero_mismatches_in_overlap
    else:
        raise ValueError('%s is not a valid merging method.' % merge_method)

    stats = FASTQMergerStats()

    # Do not use the fastqlib and fastalib objects to limit overhead.
    input1_file = gzip.open(input1_path, 'rt') if input1_is_gzipped else open(input1_path)
    input2_file = gzip.open(input2_path, 'rt') if input2_is_gzipped else open(input2_path)

    singlethreaded = True if r1_start_position == 0 and r1_end_position == -1 else False
    if r1_start_position == -1:
        return
    input1_file.seek(r1_start_position)
    input2_file.seek(r2_start_position)

    merged_file = open(merged_path, 'w')
    merge_failed_file = open(merge_failed_path, 'w')
    merge_failed_with_Ns_file = open(merge_failed_with_Ns_path, 'w')
    merge_failed_Q30_file = None if merge_failed_Q30_path is '' \
        else open(merge_failed_Q30_path, 'w')
    r1_prefix_file = open(r1_prefix_path, 'w') if r1_prefix_path else None
    r2_prefix_file = open(r2_prefix_path, 'w') if r2_prefix_path else None

    i = 0
    while True:
        if verbose:
            print_merging_progress(
                dataset_index,
                total_dataset_count,
                num_items_total,
                num_items_processed,
                next_percentage,
                num_pairs_prefix_passed,
                num_merged_pairs,
                num_merged_pairs_zero_mismatches)
            with num_items_processed.get_lock():
                num_items_processed.value += 1
        i += 1
        r1_lines = [input1_file.readline().rstrip() for _ in range(4)]
        r2_lines = [input2_file.readline().rstrip() for _ in range(4)]

        if r1_lines[0] == r1_end_string:
            # Defline strings should be unique, but double check the chunk ending position.
            if input1_file.tell() >= r1_end_position or r1_end_string == '':
                break

        stats.actual_number_of_pairs += 1

        # Check for prefix sequences.
        failed_prefix = False
        if r1_prefix_compiled:
            r1_prefix_match = (
                r1_prefix_compiled.search(r1_lines[1]) if r1_prefix_compiled else None)
            if not r1_prefix_match:
                stats.prefix_failed_in_pair_1_total += 1
                failed_prefix = True
        if r2_prefix_compiled:
            r2_prefix_match = (
                r2_prefix_compiled.search(r2_lines[1]) if r2_prefix_compiled else None)
            if not r2_prefix_match:
                stats.prefix_failed_in_pair_2_total += 1
                failed_prefix = True
            if r1_prefix_compiled:
                if not r1_prefix_match:
                    stats.prefix_failed_in_both_pairs_total += 1
        if failed_prefix:
            stats.failed_prefix_total += 1
            continue
        stats.passed_prefix_total += 1
        if verbose:
            with num_pairs_prefix_passed.get_lock():
                num_pairs_prefix_passed.value += 1

        r1_entry = FastQEntry(r1_lines, raw=ignore_deflines)
        r2_entry = FastQEntry(r2_lines, raw=ignore_deflines)

        if r1_prefix_compiled:
            r1_entry.trim(trim_from=r1_prefix_match.end())
            r1_suffix_length = read1_length - len(r1_prefix_match.group(0))
        else:
            r1_suffix_length = read1_length
        if r2_prefix_compiled:
            r2_entry.trim(trim_from=r2_prefix_match.end())
            r2_suffix_length = read2_length - len(r2_prefix_match.group(0))
        else:
            r2_suffix_length = read2_length
        min_seq_len = min(r1_suffix_length, r2_suffix_length)
        r1_suffix = r1_entry.sequence
        r2_suffix = r2_entry.sequence

        while True:
            # First, check if equal length suffixes are the same.
            if r1_suffix_length == r2_suffix_length:
                if r1_suffix == r2_suffix:
                    begin_seq = ''
                    overlap_seq = r1_suffix
                    end_seq = ''
                    merged_seq = r1_suffix
                    num_mismatches = 0
                    len_overlap = len(overlap_seq)
                    p = 0
                    recovery_dict = {'none': 0, 'r1': 0, 'r2': 0}
                    merging_done_on_complete_overlap = False
                    break

            if not partial_overlap_only:
                # Second, check if unequal length suffixes match.
                # This means that the read with the longer suffix
                # goes into the prefix and possibly the adapter of the read with the shorter suffix.
                if r1_suffix_length < r2_suffix_length:
                    rc_r2_suffix = reverse_complement(r2_suffix)
                    try:
                        rc_r2_start = rc_r2_suffix.index(r1_suffix)
                        begin_seq = rc_r2_suffix[: rc_r2_start]
                        overlap_seq = r1_suffix
                        end_seq = rc_r2_suffix[rc_r2_start + r1_suffix_length: ]
                        if retain_overlap_only:
                            merged_seq = overlap_seq
                        elif skip_suffix_trimming:
                            merged_seq = rc_r2_suffix
                        else:
                            merged_seq = overlap_seq + end_seq
                        num_mismatches = 0
                        len_overlap = len(overlap_seq)
                        p = 0
                        recovery_dict = {'none': 0, 'r1': 0, 'r2': 0}
                        merging_done_on_complete_overlap = True
                        stats.merging_done_on_complete_overlap += 1
                        break
                    except ValueError:
                        pass
                elif r2_suffix_length < r1_suffix_length:
                    rc_r2_suffix = reverse_complement(r2_suffix)
                    try:
                        r1_start = r1_suffix.index(rc_r2_suffix)
                        begin_seq = r1_suffix[: r1_start]
                        overlap_seq = rc_r2_suffix
                        end_seq = r1_suffix[r1_start + r2_suffix_length: ]
                        if retain_overlap_only:
                            merged_seq = overlap_seq
                        elif skip_suffix_trimming:
                            merged_seq = r1_suffix
                        else:
                            merged_seq = begin_seq + overlap_seq
                        num_mismatches = 0
                        len_overlap = len(overlap_seq)
                        p = 0
                        recovery_dict = {'none': 0, 'r1': 0, 'r2': 0}
                        merging_done_on_complete_overlap = True
                        stats.merging_done_on_complete_overlap += 1
                        break
                    except ValueError:
                        pass

            # Third, check for partial overlap of a long insert.
            (begin_seq, overlap_seq, end_seq), num_mismatches, recovery_dict = merge_reads(
                partial_overlap_merge_function,
                r1_entry,
                r2_entry,
                min_seq_len - 1,
                min_overlap_size=min_overlap_size)
            merged_seq = begin_seq + overlap_seq + end_seq if overlap_seq else ''
            merging_done_on_complete_overlap = False

            # find out about 'p'
            len_overlap = len(overlap_seq)
            if len_overlap == 0:
                p = 1
            else:
                p = num_mismatches / len_overlap

            if not partial_overlap_only:
                # Fourth, check for "more than complete" (really partial) overlap of a short insert.
                ((alt_begin_seq, alt_overlap_seq, alt_end_seq),
                 alt_num_mismatches,
                 alt_recovery_dict) = merge_reads(
                    partial_overlap_merge_function,
                    r1_entry,
                    r2_entry,
                    min_seq_len - 1,
                    min_overlap_size=min_overlap_size,
                    complete_overlap=True)

                if retain_overlap_only or not skip_suffix_trimming:
                    alt_merged_seq = alt_overlap_seq
                else:
                    alt_merged_seq = alt_begin_seq + alt_overlap_seq + alt_end_seq

                # find out about 'p'
                alt_len_overlap = len(alt_overlap_seq)
                if alt_len_overlap == 0:
                    alt_p = 1
                else:
                    alt_p = alt_num_mismatches / alt_len_overlap

                if alt_p < p:
                    begin_seq = alt_begin_seq
                    overlap_seq = alt_overlap_seq
                    end_seq = alt_end_seq
                    num_mismatches = alt_num_mismatches
                    recovery_dict = alt_recovery_dict
                    merged_seq = alt_merged_seq
                    len_overlap = alt_len_overlap
                    p = alt_p
                    merging_done_on_complete_overlap = True
                    stats.merging_done_on_complete_overlap += 1
            break

        if enforce_Q30_check and not merging_done_on_complete_overlap:
            if not r1_entry.Q_list:
                r1_entry.process_Q_list()
                r2_entry.process_Q_list()

            r1_passed_Q30, r1_Q30 = passes_minoche_Q30(r1_entry.Q_list[: -len_overlap])
            r2_passed_Q30, r2_Q30 = passes_minoche_Q30(r2_entry.Q_list[: -len_overlap])

            pair_passed_Q30 = r1_passed_Q30 and r2_passed_Q30

            if not r1_passed_Q30 and r2_passed_Q30:
                stats.num_Q30_fails_in_read_1 += 1
            elif r1_passed_Q30 and not r2_passed_Q30:
                stats.num_Q30_fails_in_read_2 += 1
            elif not r1_passed_Q30 and not r2_passed_Q30:
                stats.num_Q30_fails_in_both += 1

        # Generate the header line.
        header_line = (
            '%s|%o:d|m/o:%f|MR:%s|Q30:%s|CO:%d|mismatches:%d' % (
                r1_entry.header_line,
                len_overlap,
                p,
                'n=%d;r1=%d;r2=%d' % (
                    recovery_dict['none'], recovery_dict['r1'], recovery_dict['r2']),
                'n/a' if not (enforce_Q30_check and not merging_done_on_complete_overlap) \
                    else '%s=%d;%s=%d' % (
                        'p' if r1_passed_Q30 else 'f',
                        r1_Q30,
                        'p' if r2_passed_Q30 else 'f',
                        r2_Q30),
                1 if merging_done_on_complete_overlap else 0,
                num_mismatches))
        header_line = (
            '%s|%s' % (
                project_name,
                header_line.replace('_', '-')))

        # FAIL CASE ~ max_num_mismatches
        if ((max_num_mismatches >= 0 and num_mismatches > max_num_mismatches)
            or (max_num_mismatches == 0 and p == 1)):
            merge_failed_file.write('>%s\n' % header_line)
            merge_failed_file.write('%s\n' % merged_seq)
            stats.record_num_mismatches(num_mismatches, 'merge failed due to max num mismatches')
            stats.merge_failed_total += 1
            stats.pairs_eliminated_due_to_max_num_mismatches += 1
            continue

        # FAIL CASE ~ m/o
        if p > max_p or len_overlap < min_overlap_size:
            merge_failed_file.write('>%s\n' % header_line)
            merge_failed_file.write('%s\n' % merged_seq)
            stats.record_num_mismatches(num_mismatches, 'merge failed due to P value')
            stats.merge_failed_total += 1
            stats.pairs_eliminated_due_to_P += 1
            continue

        # FAIL CASE ~ N
        if not ignore_Ns:
            if 'N' in merged_seq or 'n' in merged_seq:
                merge_failed_with_Ns_file.write('>%s\n' % header_line)
                merge_failed_with_Ns_file.write('%s\n' % merged_seq)
                stats.record_num_mismatches(num_mismatches, 'merge failed due to N')
                stats.merge_failed_total += 1
                stats.pairs_eliminated_due_to_Ns += 1
                continue

        # FAIL CASE ~ Q30
        if enforce_Q30_check and not merging_done_on_complete_overlap:
            if not pair_passed_Q30:
                merge_failed_Q30_file.write('>%s\n' % header_line)
                merge_failed_Q30_file.write('%s\n' % merged_seq)
                stats.record_num_mismatches(num_mismatches, 'merge failed due to Q30')
                stats.merge_failed_total += 1
                stats.pairs_eliminated_due_to_Q30 += 1
                continue

        if verbose:
            with num_merged_pairs.get_lock():
                num_merged_pairs.value += 1
            if num_mismatches == 0:
                with num_merged_pairs_zero_mismatches.get_lock():
                    num_merged_pairs_zero_mismatches.value += 1

        # FINALLY
        merged_file.write('>%s\n' % header_line)

        if retain_overlap_only:
            merged_file.write('%s\n' % overlap_seq)
        else:
            merged_file.write('%s\n' % merged_seq)
        stats.record_num_mismatches(num_mismatches, 'merge passed')
        stats.merge_passed_total += 1

        # Record the info for the successfuly merged pair in the recovery dict.
        stats.process_recovery_dict(recovery_dict)

        # Report prefix sequences.
        if r1_prefix_file:
            r1_prefix_file.write('>%s\n' % header_line)
            r1_prefix_file.write('%s\n' % r1_prefix_match.group(0))
        if r2_prefix_file:
            r2_prefix_file.write('>%s\n' % header_line)
            r2_prefix_file.write('%s\n' % r2_prefix_match.group(0))

    input1_file.close()
    input2_file.close()
    merged_file.close()
    merge_failed_file.close()
    merge_failed_with_Ns_file.close()
    if merge_failed_Q30_file:
        merge_failed_Q30_file.close()
    if r1_prefix_file:
        r1_prefix_file.close()
    if r2_prefix_file:
        r2_prefix_file.close()

    return stats


def print_merging_progress(
    dataset_index,
    total_dataset_count,
    num_pairs_total,
    num_pairs,
    next_percentage,
    passed_prefix_total,
    num_merged_pairs_passed,
    num_zero_mismatches):
    next_percentage.acquire()
    if num_pairs.value * 100 / num_pairs_total >= next_percentage.value:
        next_percentage.value += 1
        next_percentage.release()
        percent_pairs = num_pairs.value * 100 // num_pairs_total
        try:
            sys.stderr.write(
                '\r[Merging %d of %d] %.2d%% -- (num pairs processed: %s) POK: %.1f%% :: ZM: %.1f%%'
                % (dataset_index,
                total_dataset_count,
                percent_pairs,
                big_number_pretty_print(num_pairs.value),
                passed_prefix_total.value * 100 / num_pairs.value,
                num_zero_mismatches.value * 100 / num_merged_pairs_passed.value))
        except ZeroDivisionError:
            pass
        sys.stderr.flush()
    else:
        next_percentage.release()
    return


def merge_reads(
    merge_function,
    r1_entry,
    r2_entry,
    max_overlap_size,
    min_overlap_size=16,
    complete_overlap=False,
    min_qual_score=15):

    if complete_overlap:
        seq1 = reverse_complement(r2_entry.sequence)
        seq2 = r1_entry.sequence
    else:
        seq1 = r1_entry.sequence
        seq2 = reverse_complement(r2_entry.sequence)

    seq1_overlap_start, seq2_overlap_end = merge_function(
        seq1=seq1, seq2=seq2, max_overlap_size=max_overlap_size, min_overlap_size=min_overlap_size)

    begin_seq = seq1[: seq1_overlap_start].lower()
    overlap_seq1 = seq1[seq1_overlap_start: ].lower()
    overlap_seq2 = seq2[: seq2_overlap_end].lower()
    end_seq = seq2[seq2_overlap_end: ].lower()
    overlap_seq = ''

    recovery_dict = {'none': 0, 'r1': 0, 'r2': 0}
    for i in range(len(overlap_seq1)):
        if overlap_seq1[i] != overlap_seq2[i]:
            # MISMATCH FOUND!
            if complete_overlap:
                read_id_with_better_base_qual = get_read_id_with_better_base_qual(
                    r1_entry, r2_entry, i, seq1_overlap_start + i, min_qual_score)

                if read_id_with_better_base_qual == 0:
                    overlap_seq += 'N'
                    recovery_dict['none'] += 1
                elif read_id_with_better_base_qual == 1:
                    overlap_seq += overlap_seq2[i].upper()
                    recovery_dict['r1'] += 1
                elif read_id_with_better_base_qual == 2:
                    overlap_seq += overlap_seq1[i].upper()
                    recovery_dict['r2'] += 1
            else:
                read_id_with_better_base_qual = get_read_id_with_better_base_qual(
                    r1_entry, r2_entry, seq1_overlap_start + i, i, min_qual_score)

                if read_id_with_better_base_qual == 0:
                    overlap_seq += 'N'
                    recovery_dict['none'] += 1
                elif read_id_with_better_base_qual == 1:
                    overlap_seq += overlap_seq1[i].upper()
                    recovery_dict['r1'] += 1
                elif read_id_with_better_base_qual == 2:
                    overlap_seq += overlap_seq2[i].upper()
                    recovery_dict['r2'] += 1
        else:
            overlap_seq += overlap_seq1[i]

    return (begin_seq, overlap_seq, end_seq), sum(recovery_dict.values()), recovery_dict


def merge_with_zero_mismatches_in_overlap(
    seq1=None, seq2=None, max_overlap_size=None, min_overlap_size=16):
    # Go from most to least overlap, returning the first overlap found without mismatches.
    for overlap_size in range(max_overlap_size, min_overlap_size - 1, -1):
        if seq1[-overlap_size: ] == seq2[: overlap_size]:
            seq1_overlap_start = len(seq1) - overlap_size
            seq2_overlap_end = overlap_size
            break
    else:
        seq1_overlap_start = len(seq1)
        seq2_overlap_end = 0
    return seq1_overlap_start, seq2_overlap_end



def merge_by_distance_metric(
    metric, seq1=None, seq2=None, max_overlap_size=None, min_overlap_size=16):
    smallest_dist = sys.maxsize
    selected_overlap_size = 0
    # Go from least to most overlap, since the entire overlap range is tested.
    for overlap_size in range(min_overlap_size, max_overlap_size + 1):
        dist = metric(seq1[-overlap_size: ], seq2[: overlap_size])
        if dist <= smallest_dist:
            smallest_dist = dist
            selected_overlap_size = overlap_size
    seq1_overlap_start = len(seq1) - selected_overlap_size
    seq2_overlap_end = selected_overlap_size
    return seq1_overlap_start, seq2_overlap_end


def get_read_id_with_better_base_qual(
    r1_entry, r2_entry, base_index_in_read_1, base_index_in_reversed_read_2, min_qual_score):
    """
    When there is a disagreement regarding a base between pairs at the overlapped region,
    such as this one:

    read 1: AAAAAAAAAAAAAAAAAAAAAAAAAAAA
    read 2:                AAAATAAAAAAAAAAAAAAAAAAAAAAAA
                               ^

    there can be two outcomes for the merged sequence:
    We would either use the base in read 1,
    in which case the merged sequence would look like this:

    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                       ^

    or we would use the base in read 2 that would result with this one:

    AAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAA
                       ^


    This function takes both reads and returns the read id (1 or 2)
    that has a higher quality score assigned by the sequencer for the base in question.
    If the higher quality score is still lower than min_qual_score parameter,
    0 is returned as read id.
    """

    if not r1_entry.Q_list:
        r1_entry.process_Q_list()
        r2_entry.process_Q_list()

    base_qual_in_read_1 = r1_entry.Q_list[base_index_in_read_1]
    base_qual_in_read_2 = r2_entry.Q_list[:: -1][base_index_in_reversed_read_2]

    # If neither of the bases in question satisfies min_qual_score expectation, return 0.
    # In this case, the dispute resolves to an ambiguous base in the merged sequence.
    if max([base_qual_in_read_1, base_qual_in_read_2]) < min_qual_score:
        return 0

    if base_qual_in_read_1 >= base_qual_in_read_2:
        return 1
    else:
        return 2


def passes_minoche_Q30(base_qualities):
    # This algorithm is from Minoche et al.
    # It calculates the length of the read,
    # returning True only if two-thirds of bases
    # in the first half of the read have Q-scores over Q30.

    half_length = len(base_qualities) / 2
    Q30 = len([True for _q in base_qualities[: int(half_length)] if _q > 30])
    if Q30 < (0.66 * half_length):
        return (False, Q30)
    return True, Q30
