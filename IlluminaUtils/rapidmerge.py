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

import gzip
import multiprocessing
import os
import re
import stat
import struct

from datetime import datetime

from IlluminaUtils.utils.helperfunctions import (
    combine_files, conv_dict, is_file_exists, reverse_complement)
from IlluminaUtils.utils.terminal import Progress
from IlluminaUtils.lib.fastqlib import FastQSource, FastQEntry


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
        input1_is_gzipped=False,
        input2_is_gzipped=False,
        ignore_deflines=False,
        output_dir='',
        output_file_name='output',
        r1_prefix_pattern='',
        r2_prefix_pattern='',
        report_r1_prefix=False,
        report_r2_prefix=False,
        min_overlap_size=16,
        partial_overlap_only=True,
        retain_overlap_only=False,
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
        assert len(first_seq1) == len(first_seq2)
        self.min_seq_len = min(
            len(first_seq1) - len(r1_prefix_pattern), len(first_seq2) - len(r1_prefix_pattern))

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
                self.min_seq_len,
                input1_is_gzipped=self.input1_is_gzipped,
                input2_is_gzipped=self.input2_is_gzipped,
                ignore_deflines=self.ignore_deflines,
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
        # Find positions at which to chunk the input files.
        start_positions, end_positions, end_strings = self.find_fastq_chunk_starts(
            self.input1_path, self.num_cores, self.input1_is_gzipped)
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
                *(
                    merge_method,
                    self.input1_path,
                    self.input2_path,
                    self.min_seq_len),
                **{
                    'input1_is_gzipped': self.input1_is_gzipped,
                    'input2_is_gzipped': self.input2_is_gzipped,
                    'ignore_deflines': self.ignore_deflines,
                    'merged_path': temp_merged_path,
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


    def find_fastq_chunk_starts(self, path, num_chunks, input_is_gzipped):
        if input_is_gzipped:
            # Uncompressed size is stored in the last four digits of a gzip file.
            with open(path, 'rb') as f:
                f.seek(0, 2)
                f.seek(-4, 2)
                uncompressed_file_size = struct.unpack('I', f.read(4))[0]
        else:
            uncompressed_file_size = os.stat(path)[stat.ST_SIZE]
        chunk_size = uncompressed_file_size // num_chunks

        fastq_file = gzip.open(path, 'rt') if input_is_gzipped else open(path)
        start_positions = []
        end_strings = []
        position = 0
        prev_position = 0
        for chunk in range(num_chunks):
            fastq_file.seek(position)
            while True:
                line_end = fastq_file.readline()
                # Checking for empty string must come before checking for '@'
                if line_end == '':
                    # If EOF, append -1 rather than the last position.
                    start_positions.append(-1)
                    break
                prev_position = position
                position = fastq_file.tell()
                if line_end[0] == '@':
                    start_positions.append(prev_position)
                    position += chunk_size
                    break
            if chunk > 0:
                end_strings.append(line_end.rstrip())
        end_strings.append('')
        end_positions = [p for p in start_positions[1:]] + [-1]
        fastq_file.close()
        return start_positions, end_positions, end_strings


def merge_reads_in_files(
    merge_method,
    input1_path,
    input2_path,
    min_seq_len,
    input1_is_gzipped=False,
    input2_is_gzipped=False,
    ignore_deflines=False,
    merged_path='output_MERGED',
    r1_prefix_compiled=None,
    r2_prefix_compiled=None,
    r1_prefix_path='',
    r2_prefix_path='',
    min_overlap_size=16,
    partial_overlap_only=True,
    retain_overlap_only=False,
    start_position=0,
    end_position=-1,
    end_string=''):

    total_pairs_count = 0
    merged_count = 0
    fully_overlapping_count = 0
    prefix_passed_count = 0
    total_prefix_failed_count = 0
    r1_prefix_failed_count = 0
    r2_prefix_failed_count = 0
    both_prefixes_failed_count = 0
    pair_disqualified_by_Ns_count = 0

    # Do not use the fastqlib and fastalib objects to limit overhead.
    input1_file = gzip.open(input1_path, 'rt') if input1_is_gzipped else open(input1_path)
    input2_file = gzip.open(input2_path, 'rt') if input2_is_gzipped else open(input2_path)

    if start_position == -1:
        return
    input1_file.seek(start_position)
    input2_file.seek(start_position)

    merged_file = open(merged_path, 'w')
    r1_prefix_file = open(r1_prefix_path, 'w') if r1_prefix_path else None
    r2_prefix_file = open(r2_prefix_path, 'w') if r2_prefix_path else None

    while True:
        r1_lines = [input1_file.readline().rstrip() for _ in range(4)]
        r2_lines = [input2_file.readline().rstrip() for _ in range(4)]

        if r1_lines[0] == end_string:
            # Defline strings should be unique, but double check the chunk ending position.
            if input1_file.tell() >= end_position or end_string == '':
                break

        total_pairs_count += 1

        # Check for prefix sequences.
        if r1_prefix_compiled:
            r1_prefix_match = r1_prefix_compiled.search(r1_lines[1]) if r1_prefix_compiled else None
            if not r1_prefix_match:
                r1_prefix_failed_count += 1
        if r2_prefix_compiled:
            r2_prefix_match = r2_prefix_compiled.search(r2_lines[1]) if r2_prefix_compiled else None
            if not r2_prefix_match:
                r2_prefix_failed_count += 1
            if r1_prefix_compiled:
                if not r1_prefix_match:
                    both_prefixes_failed_count += 1
        if r1_prefix_compiled:
            if not r1_prefix_match:
                total_prefix_failed_count += 1
                continue
        if r2_prefix_compiled:
            if not r2_prefix_match:
                total_prefix_failed_count += 1
                continue
        prefix_passed_count += 1

        r1_entry = FastQEntry(r1_lines, raw=ignore_deflines)
        r2_entry = FastQEntry(r2_lines, raw=ignore_deflines)

        if r1_prefix_compiled:
            r1_entry.trim(trim_from=r1_prefix_match.end())
        if r2_prefix_compiled:
            r2_entry.trim(trim_from=r2_prefix_match.end())

        begin_seq, overlap_seq, end_seq = merge_method(
            r1_entry.sequence,
            r2_entry.sequence,
            min_seq_len,
            min_overlap_size=min_overlap_size)
        merged_seq = begin_seq + overlap_seq + end_seq if overlap_seq else ''
        partially_overlapping = True
        if not partial_overlap_only:
            alt_begin_seq, alt_overlap_seq, alt_end_seq = merge_method(
                r1_entry.sequence,
                r2_entry.sequence,
                min_seq_len,
                min_overlap_size=min_overlap_size,
                complete_overlap=True)
            if len(alt_overlap_seq) > len(merged_seq):
                merged_seq = alt_overlap_seq
                begin_seq = alt_begin_seq
                overlap_seq = alt_overlap_seq
                end_seq = alt_end_seq
        overlap_size = len(overlap_seq)

        insert_seq = merged_seq

        if not insert_seq:
            continue
        if 'N' in insert_seq:
            pair_disqualified_by_Ns_count += 1
            continue
        merged_count += 1
        if not partially_overlapping:
            fully_overlapping_count += 1

        iu_formatted_defline = (">%s|o:%d|m/o:0|MR:n=0;r1=0;r2=0|Q30:n/a|mismatches:0"
                                % (r1_entry.header_line, overlap_size))
        merged_file.write("%s\n%s\n" % (iu_formatted_defline, insert_seq))

        # Report prefix sequences.
        if r1_prefix_file:
            r1_prefix_file.write('>' + iu_formatted_defline)
            r1_prefix_file.write(r1_prefix_match.group(0) + '\n')
        if r2_prefix_file:
            r2_prefix_file.write('>' + iu_formatted_defline)
            r2_prefix_file.write(r2_prefix_match.group(0) + '\n')

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
        'total_prefix_failed': total_prefix_failed_count,
        'r1_prefix_failed': r1_prefix_failed_count,
        'r2_prefix_failed': r2_prefix_failed_count,
        'both_prefixes_failed': both_prefixes_failed_count,
        'pair_disqualified_by_Ns': pair_disqualified_by_Ns_count}

    return count_stats


def merge_with_zero_mismatches_in_overlap(
    r1_seq,
    r2_seq,
    max_overlap_size,
    min_overlap_size=16,
    complete_overlap=False):

    if complete_overlap:
        seq1 = reverse_complement(r2_seq)
        seq2 = r1_seq
    else:
        seq1 = r1_seq
        seq2 = reverse_complement(r2_seq)

    for overlap_size in range(max_overlap_size, min_overlap_size - 1, -1):
        if seq1[-overlap_size: ] == seq2[: overlap_size]:
            overlap_start1 = len(seq1) - overlap_size
            overlap_end2 = overlap_size
            break
    else:
        overlap_start1 = len(seq1)
        overlap_end2 = 0
    begin_seq = seq1[: overlap_start1]
    overlap_seq = seq1[overlap_start1: ]
    end_seq = seq2[overlap_end2: ]

    return begin_seq, overlap_seq, end_seq