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

import multiprocessing
import os
import re
import stat

from datetime import datetime

from IlluminaUtils.utils.helperfunctions import (
    combine_files, conv_dict, is_file_exists, reverse_complement)
from IlluminaUtils.utils.terminal import Progress


class Multiprocessor:
    def __init__(self, num_cores, task_name=''):
        self.num_cores = num_cores
        self.task_name = task_name
        self.pool = multiprocessing.Pool(num_cores)
        self.jobs = []
        return


    def start_job(self, process, *args, **kwargs):
        self.jobs.append(self.pool.apply_async(process, args, kwargs))
        return


    def get_output(self, verbose=True):
        outputs = []
        completed_jobs = []
        if verbose:
            progress = Progress()
            progress.new(os.getpid())
            progress.update(
                "%d of %d %s jobs complete" % (len(outputs), len(self.jobs), self.task_name))
        while True:
            for job in self.jobs:
                if job in completed_jobs:
                    continue
                if job.ready():
                    outputs.append(job.get())
                    completed_jobs.append(job)
                    if verbose:
                        progress.update("%d of %d %s jobs complete"
                                    % (len(completed_jobs), len(self.jobs), self.task_name))
            if len(completed_jobs) == len(self.jobs):
                break
        self.pool.close()
        if verbose:
            print()
        return outputs


class FastQMerger:
    def __init__(
        self,
        input1_path,
        input2_path,
        output_dir='',
        output_file_name='output',
        r1_prefix_pattern='',
        r2_prefix_pattern='',
        report_r1_prefix=False,
        report_r2_prefix=False,
        min_overlap_size=15,
        partial_overlap_only=True,
        retain_overlap_only=False,
        num_cores=1):

        is_file_exists(input1_path)
        is_file_exists(input2_path)
        self.input1_path = input1_path
        self.input2_path = input2_path

        output_path_maker = lambda p: os.path.join(output_dir, output_file_name + p)
        self.merged_path = output_path_maker('_MERGED')

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

        self.min_overlap_size = min_overlap_size
        self.partial_overlap_only = partial_overlap_only
        self.retain_overlap_only = retain_overlap_only

        if not 0 < num_cores <= multiprocessing.cpu_count():
            raise RuntimeError("\"%d\" is not a valid number of cores. "
                               "The number of cores be between 1 and %d."
                               % (num_cores, multiprocessing.cpu_count()))
        self.num_cores = num_cores
        return


    def run(self, merge_method):
        # singlethreaded
        if self.num_cores == 1:
            count_stats = merge_reads_in_files(
                merge_method,
                self.input1_path,
                self.input2_path,
                merged_path=self.merged_path,
                r1_prefix_compiled=self.r1_prefix_compiled,
                r2_prefix_compiled=self.r2_prefix_compiled,
                r1_prefix_path=self.r1_prefix_path,
                r2_prefix_path=self.r2_prefix_path,
                min_overlap_size=self.min_overlap_size,
                partial_overlap_only=self.partial_overlap_only,
                retain_overlap_only=self.retain_overlap_only)
            return count_stats

        # Prepare multiprocessing.
        progress = Progress()
        progress.new(os.getpid())
        progress.update("Setting up read merging jobs")
        print()
        # Break input files into chunks of roughly the same numbers of lines.
        start_positions, end_strings = self.find_fastq_chunk_starts(
            self.input1_path, self.num_cores)
        end_positions = [next_start_position for next_start_position in start_positions[1:]] + [-1]
        # Create temporary output files.
        time_str = datetime.now().isoformat(timespec='seconds').replace('-', '').replace(':', '')
        temp_merged_paths = [
            ("%s_TEMP_%d_%s"
             % (self.merged_path, chunk_index, time_str))
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

        # Spawn jobs.
        multiprocessor = Multiprocessor(self.num_cores, task_name='read merging')
        for (temp_merged_path,
             temp_r1_prefix_path,
             temp_r2_prefix_path,
             start_position,
             end_position,
             end_string) in zip(temp_merged_paths,
                                temp_r1_prefix_paths,
                                temp_r2_prefix_paths,
                                start_positions,
                                end_positions,
                                end_strings):
            multiprocessor.start_job(
                merge_reads_in_files,
                *(merge_method,
                 self.input1_path,
                 self.input2_path),
                **{'merged_path': temp_merged_path,
                 'r1_prefix_compiled': self.r1_prefix_compiled,
                 'r2_prefix_compiled': self.r2_prefix_compiled,
                 'r1_prefix_path': temp_r1_prefix_path,
                 'r2_prefix_path': temp_r2_prefix_path,
                 'min_overlap_size': self.min_overlap_size,
                 'partial_overlap_only': self.partial_overlap_only,
                 'retain_overlap_only': self.retain_overlap_only,
                 'start_position': start_position,
                 'end_position': end_position,
                 'end_string': end_string})
        count_stats_chunks = multiprocessor.get_output()

        # Delete temp files after combining them.
        progress.update("Combining temporary files produced by each job")
        print()
        progress.end()
        combine_files(temp_merged_paths, self.merged_path)
        [os.remove(temp_merged_path) for temp_merged_path in temp_merged_paths]
        if self.r1_prefix_path:
            combine_files(temp_r1_prefix_paths, self.r1_prefix_path)
            [os.remove(temp_r1_prefix_path) for temp_r1_prefix_path in temp_r1_prefix_paths]
        if self.r2_prefix_path:
            combine_files(temp_r2_prefix_paths, self.r2_prefix_path)
            [os.remove(temp_r2_prefix_path) for temp_r2_prefix_path in temp_r2_prefix_paths]

        # Sum counts from each chunk.
        count_stats = dict([(key, 0) for key in count_stats_chunks[0]])
        for count_stats_chunk in count_stats_chunks:
            for key in count_stats_chunk:
                count_stats[key] += count_stats_chunk[key]
        return count_stats


    def find_fastq_chunk_starts(self, path, num_chunks):
        fastq_file = open(path)
        file_size = os.stat(path)[stat.ST_SIZE]
        chunk_size = file_size // num_chunks
        chunk_start_positions = []
        chunk_end_strings = []
        position = 0
        prev_position = 0
        for chunk in range(num_chunks):
            fastq_file.seek(position)
            while True:
                line_end = fastq_file.readline()
                # Checking for empty string must come before checking for '@'
                if line_end == '':
                    # If EOF, append -1 rather than the last position.
                    chunk_start_positions.append(-1)
                    break
                prev_position = position
                position = fastq_file.tell()
                if line_end[0] == '@':
                    chunk_start_positions.append(prev_position)
                    position += chunk_size
                    break
            if chunk > 0:
                chunk_end_strings.append(line_end)
        chunk_end_strings.append('')
        fastq_file.close()
        return chunk_start_positions, chunk_end_strings


def merge_reads_in_files(
    merge_method,
    input1_path,
    input2_path,
    merged_path='output_MERGED',
    r1_prefix_compiled=None,
    r2_prefix_compiled=None,
    r1_prefix_path='',
    r2_prefix_path='',
    min_overlap_size=15,
    partial_overlap_only=True,
    retain_overlap_only=False,
    start_position=0,
    end_position=-1,
    end_string=''):

    total_pairs_count = 0
    merged_count = 0
    fully_overlapping_count = 0
    prefix_passed_count = 0
    r1_prefix_failed_count = 0
    r2_prefix_failed_count = 0
    both_prefixes_failed_count = 0
    pair_disqualified_by_Ns_count = 0

    # Do not use the fastqlib and fastalib objects to limit overhead.
    input1_file = open(input1_path)
    input2_file = open(input2_path)
    if start_position == -1:
        return
    input1_file.seek(start_position)
    input2_file.seek(start_position)

    merged_file = open(merged_path, 'w')
    r1_prefix_file = open(r1_prefix_path, 'w') if r1_prefix_path else None
    r2_prefix_file = open(r2_prefix_path, 'w') if r2_prefix_path else None

    block_index = 0
    while True:
        r1_line = input1_file.readline()
        r2_line = input2_file.readline()

        if r1_line == end_string:
            if input1_file.tell() >= end_position or end_string == '':
                break

        if block_index == 0:
            r1_defline = r1_line
        elif block_index == 1:
            r1_seq = r1_line.rstrip()
            r2_seq = r2_line.rstrip()
        elif block_index == 3:
            # Reset the sequence block.
            block_index = 0
            total_pairs_count += 1

            # Check for prefix sequences.
            r1_prefix_match = r1_prefix_compiled.search(r1_seq) if r1_prefix_compiled else None
            r2_prefix_match = r2_prefix_compiled.search(r2_seq) if r2_prefix_compiled else None
            if not r1_prefix_match:
                r1_prefix_failed_count += 1
            if not r2_prefix_match:
                r2_prefix_failed_count += 1
                if not r1_prefix_match:
                    both_prefixes_failed_count += 1
            if not (r1_prefix_match and r2_prefix_match):
                continue
            else:
                prefix_passed_count += 1
            r1_prefix_seq = r1_prefix_match.group(0) if r1_prefix_match else ''
            r2_prefix_seq = r2_prefix_match.group(0) if r2_prefix_match else ''

            insert_seq, overlap_size, partially_overlapping = merge_method(
                r1_defline,
                r1_seq,
                r2_seq,
                r1_prefix_length=len(r1_prefix_seq),
                r2_prefix_length=len(r2_prefix_seq),
                min_overlap_size=min_overlap_size,
                partial_overlap_only=partial_overlap_only,
                retain_overlap_only=retain_overlap_only)

            if not insert_seq:
                continue
            if 'N' in insert_seq:
                pair_disqualified_by_Ns_count += 1
                continue
            merged_count += 1
            if not partially_overlapping:
                fully_overlapping_count += 1

            iu_formatted_defline = (">%s|o:%d|m/o:0|MR:n=0;r1=0;r2=0|Q30:n/a|mismatches:0"
                                    % (r1_defline.rstrip()[1:], overlap_size))
            merged_file.write("%s\n%s\n" % (iu_formatted_defline, insert_seq))

            # Report prefix sequences.
            if r1_prefix_file:
                r1_prefix_file.write('>' + iu_formatted_defline)
                r1_prefix_file.write(r1_prefix_seq + '\n')
            if r2_prefix_file:
                r2_prefix_file.write('>' + iu_formatted_defline)
                r2_prefix_file.write(r2_prefix_seq + '\n')
            continue

        # This increment is not reached after the last line in the sequence block.
        block_index += 1

    input1_file.close()
    input2_file.close()
    merged_file.close()
    if r1_prefix_file:
        r1_prefix_file.close()
    if r2_prefix_file:
        r2_prefix_file.close()

    count_stats = {
        'total_pairs': total_pairs_count,
        'merged': merged_count,
        'fully_overlapping_count': fully_overlapping_count,
        'prefix_passed': prefix_passed_count,
        'r1_prefix_failed': r1_prefix_failed_count,
        'r2_prefix_failed': r2_prefix_failed_count,
        'both_prefixes_failed': both_prefixes_failed_count,
        'pair_disqualified_by_Ns': pair_disqualified_by_Ns_count}

    return count_stats


def merge_with_zero_mismatches_in_overlap(
    r1_defline,
    r1_seq,
    r2_seq,
    r1_prefix_length=0,
    r2_prefix_length=0,
    min_overlap_size=15,
    partial_overlap_only=True,
    retain_overlap_only=False):
    # Returns insert sequence without prefix or suffix sequences, insert length,
    # and Boolean indicating partial (True), full (False), or no (None) overlap.

    # Reads 1 and 2 are assumed to be of equal length.
    assert len(r1_seq) == len(r2_seq)
    read_length = len(r1_seq)
    if read_length <= min_overlap_size:
        return '', 0, None

    r2_rc_seq = reverse_complement(r2_seq)

    if not partial_overlap_only:
        full_overlap_r2_end_index = read_length - r2_prefix_length

        # Check if the insert is the same size as the read.
        r1_end_index = read_length - r2_prefix_length
        r2_start_index = r1_prefix_length
        if (r1_seq[r1_prefix_length: r1_end_index]
            == r2_rc_seq[r2_start_index: full_overlap_r2_end_index]):
            return r1_seq[r1_prefix_length: r1_end_index], r1_end_index - r1_prefix_length, False

    shift = 1
    while read_length - shift > min_overlap_size:
        # Check for partial overlap.
        if read_length - shift > min_overlap_size:
            if shift > r1_prefix_length:
                r1_start_index = shift
                r2_start_index = 0
            else:
                r1_start_index = r1_prefix_length
                r2_start_index = r1_prefix_length - shift
            r1_end_index = (read_length - r2_prefix_length + shift if shift - r2_prefix_length < 0
                            else read_length)
            r2_end_index = (read_length - shift if shift > r2_prefix_length
                            else read_length - r2_prefix_length)
            if (r1_seq[r1_start_index: r1_end_index]
                == r2_rc_seq[r2_start_index: r2_end_index]):
                if retain_overlap_only:
                    insert_seq = r1_seq[r1_start_index: r1_end_index]
                else:
                    insert_seq = (r1_seq[r1_prefix_length: r1_end_index]
                                  + r2_rc_seq[r2_end_index: read_length - r2_prefix_length])
                return insert_seq, r1_end_index - r1_start_index, True

        # Check for full overlap.
        if partial_overlap_only:
            shift += 1
            continue
        if read_length - shift - r1_prefix_length - r2_prefix_length <= min_overlap_size:
            shift += 1
            continue
        r1_end_index = read_length - r2_prefix_length - shift
        r2_start_index = r1_prefix_length + shift
        if (r1_seq[r1_prefix_length: r1_end_index]
            == r2_rc_seq[r2_start_index: full_overlap_r2_end_index]):
            return r1_seq[r1_prefix_length: r1_end_index], r1_end_index - r1_prefix_length, False
        shift += 1
    return '', 0, None


### PYTEST TESTS
def test_merge_read1_read2_fully_overlapping():
    r1_seq = 'AAAAACCCCCGGGGGTTTTT'
    r2_seq = 'CCCCCGGGGGTTTTTAAAAA' # rc = TTTTTAAAAACCCCCGGGGG

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_prefix_length=0,
        r2_prefix_length=0,
        partial_overlap_only=False) == (
            'AAAAACCCCCGGGGG', 'EEEEEEEEEEEEEEE', 15, False)


def test_merge_read1_read2_partially_overlapping():
    r1_seq = 'AAAAATTTTTGGGGGCCCCC'
    r2_seq = 'TTTTTGGGGGCCCCCAAAAA' # rc = TTTTTGGGGGCCCCCAAAAA

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_prefix_length=0,
        r2_prefix_length=0,
        partial_overlap_only=False) == (
            'AAAAATTTTTGGGGGCCCCCAAAAA', 'JJJJJEEEEEEEEEEEEEEEAAAAA', 15, True)