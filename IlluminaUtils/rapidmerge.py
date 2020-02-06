#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
import os
import re
import stat

from datetime import datetime

from IlluminaUtils.utils.helperfunctions import (
    combine_files, conv_dict, is_file_exists, reverse_complement)
from IlluminaUtils.utils.terminal import Progress


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
        check_full_overlap=False,
        retain_only_overlap=False,
        trim_suffix=False,
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
        self.check_full_overlap = check_full_overlap
        self.retain_only_overlap = retain_only_overlap
        self.trim_suffix = trim_suffix

        if not 0 < num_cores <= multiprocessing.cpu_count():
            raise RuntimeError("\"%d\" is not a valid number of cores. "
                               "The number of cores be between 1 and %d."
                               % (num_cores, multiprocessing.cpu_count()))
        self.num_cores = num_cores
        return


    def run(self):
        # singlethreaded
        if self.num_cores == 1:
            count_stats = merge_reads_in_files(
                input1_path=self.input1_path,
                input2_path=self.input2_path,
                merged_path=self.merged_path,
                r1_prefix_compiled=self.r1_prefix_compiled,
                r2_prefix_compiled=self.r2_prefix_compiled,
                r1_prefix_path=self.r1_prefix_path,
                r2_prefix_path=self.r2_prefix_path,
                min_overlap_size=self.min_overlap_size,
                check_full_overlap=self.check_full_overlap,
                retain_only_overlap=self.retain_only_overlap,
                trim_suffix=self.trim_suffix)
            return count_stats

        # Prepare multiprocessing.
        progress = Progress()
        progress.new(os.getpid())
        progress.update("Setting up read merging jobs")
        print()
        progress.end()
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
        pool = multiprocessing.Pool(self.num_cores)
        jobs = []
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
            jobs.append(pool.apply_async(
                merge_reads_in_files,
                (self.input1_path,
                 self.input2_path),
                {'merged_path': temp_merged_path,
                 'r1_prefix_compiled': self.r1_prefix_compiled,
                 'r2_prefix_compiled': self.r2_prefix_compiled,
                 'r1_prefix_path': temp_r1_prefix_path,
                 'r2_prefix_path': temp_r2_prefix_path,
                 'min_overlap_size': self.min_overlap_size,
                 'check_full_overlap': self.check_full_overlap,
                 'retain_only_overlap': self.retain_only_overlap,
                 'trim_suffix': self.trim_suffix,
                 'start_position': start_position,
                 'end_position': end_position,
                 'end_string': end_string}))

        # Wait for jobs to finish.
        count_stats_chunks = []
        completed_jobs = []
        progress.new(os.getpid())
        progress.update(
            "%d of %d read merging jobs complete" % (len(count_stats_chunks), len(jobs)))
        while True:
            for job in jobs:
                if job not in completed_jobs:
                    if job.ready():
                        count_stats_chunks.append(job.get())
                        completed_jobs.append(job)
                        progress.update("%d of %d read merging jobs complete"
                                    % (len(completed_jobs), len(jobs)))
            if len(completed_jobs) == len(jobs):
                break
        pool.close()
        print()
        progress.end()

        # Delete temp files after combining them.
        progress.new(os.getpid())
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
    input1_path,
    input2_path,
    merged_path='output_MERGED',
    r1_prefix_compiled=None,
    r2_prefix_compiled=None,
    r1_prefix_path='',
    r2_prefix_path='',
    min_overlap_size=15,
    check_full_overlap=False,
    retain_only_overlap=False,
    trim_suffix=False,
    start_position=0,
    end_position=-1,
    end_string=''):

    total_pairs_count = 0
    merged_count = 0
    merged_r1_does_not_precede_r2_count = 0
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
            r2_defline = r2_line
        elif block_index == 1:
            r1_seq = r1_line.rstrip()
            r2_seq = r2_line.rstrip()
        elif block_index == 3:
            # Reset the sequence block.
            block_index = 0
            total_pairs_count += 1

            failure_before_merge = False
            if 'N' in r1_seq or 'N' in r2_seq:
                pair_disqualified_by_Ns_count += 1
                failure_before_merge = True

            r1_score = r1_line.rstrip()
            r2_score = r2_line.rstrip()
            # Use the read 1 defline for merged reads (arbitrary).
            defline = r1_defline

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
                failure_before_merge = True
            else:
                prefix_passed_count += 1

            if failure_before_merge:
                continue

            r1_prefix_seq = r1_prefix_match.group(0) if r1_prefix_match else ''
            r2_prefix_seq = r2_prefix_match.group(0) if r2_prefix_match else ''

            insert_seq, insert_score, overlap_size, r1_precedes_r2 = merge_read1_read2(
                r1_seq,
                r2_seq,
                r1_score,
                r2_score,
                r1_prefix_length=len(r1_prefix_seq),
                r2_prefix_length=len(r2_prefix_seq),
                min_overlap_size=min_overlap_size,
                check_full_overlap=check_full_overlap,
                retain_only_overlap=retain_only_overlap,
                trim_suffix=trim_suffix)

            if not insert_seq:
                continue
            merged_count += 1
            if not r1_precedes_r2:
                merged_r1_does_not_precede_r2_count += 1

            iu_formatted_defline = (">%s|o:%d|m/o:0|MR:n=0;r1=0;r2=0|Q30:n/a|mismatches:0"
                                    % (defline.rstrip()[1:], overlap_size))
            merged_file.write("%s\n%s\n" % (iu_formatted_defline, insert_seq))

            # Report prefix sequences.
            if r1_prefix_file:
                r1_prefix_file.write('>' + defline[1:])
                r1_prefix_file.write(r1_prefix_seq + '\n')
            if r2_prefix_file:
                r2_prefix_file.write('>' + defline[1:])
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
        'merged_r1_does_not_precede_r2': merged_r1_does_not_precede_r2_count,
        'prefix_passed': prefix_passed_count,
        'r1_prefix_failed': r1_prefix_failed_count,
        'r2_prefix_failed': r2_prefix_failed_count,
        'both_prefixes_failed': both_prefixes_failed_count,
        'pair_disqualified_by_Ns': pair_disqualified_by_Ns_count}

    return count_stats


def merge_read1_read2(
    r1_seq,
    r2_seq,
    r1_score,
    r2_score,
    r1_prefix_length=0,
    r2_prefix_length=0,
    min_overlap_size=15,
    check_full_overlap=False,
    retain_only_overlap=False,
    trim_suffix=False):

    r2_rc_seq = reverse_complement(r2_seq)
    r1_overlap_start, r2_overlap_start, overlap_size = find_longest_read_overlap(
        seq1=r1_seq,
        seq2=r2_rc_seq,
        min_overlap_size=min_overlap_size,
        seq2_can_precede_seq1=check_full_overlap)
    if overlap_size == 0:
        insert_seq = ''
        insert_score = ''
        return insert_seq, insert_score, overlap_size, None
    r1_overlap_end = r1_overlap_start + overlap_size
    r2_overlap_end = r2_overlap_start + overlap_size

    # FULL OVERLAP OF INSERT
    # r1_seq:         CCCCCGGGGGtacgt
    # r2_rc_seq: acgtaCCCCCGGGGG
    # Here, r1 reads into r2 adapter (tacgt...), and r2 reads into r1 adapter (...acgta).
    # Overlap needs to start at index 0 of r1 and end at index -1 of r2;
    # if r1 index start > 0 or r2 index end < -1, there are mismatches.

    # PARTIAL OVERLAP OF INSERT
    # r1_seq:    AAAAACCCCCGGGGG
    # r2_rc_seq:      CCCCCGGGGGTTTTT
    # Overlap needs to end at index -1 of r1 and start at index 0 of r2;
    # if r1 index end < -1 or r2 index start > 0, there are mismatches.

    r1_precedes_r2 = None
    # Overlap scores are from both read 1 and read 2; the scheme was arbitrarily chosen.
    # FULL OVERLAP
    if r1_overlap_start == 0 and r2_overlap_end == len(r2_rc_seq):
        # Only overlapping sequence is retained regardless of `retain_only_overlap`.
        r1_insert_start = r1_overlap_start + r1_prefix_length if trim_suffix else r1_overlap_start
        r1_insert_end = r1_overlap_end - r2_prefix_length if trim_suffix else r1_overlap_end
        insert_seq = r1_seq[r1_insert_start: r1_insert_end]
        insert_score = r1_score[r1_insert_start: r1_insert_end]
        r1_precedes_r2 = False
    # PARTIAL OVERLAP
    elif r1_overlap_end == len(r1_seq) and r2_overlap_start == 0:
        if retain_only_overlap:
            r1_insert_start = (max(r1_prefix_length, r1_overlap_start)
                               if trim_suffix else r1_overlap_start)
            r2_insert_end = (min(len(r2_rc_seq) - r2_prefix_length, r2_overlap_end)
                             if trim_suffix else r2_overlap_end)
        else:
            # To reconstruct the insert,
            # (arbitrarily) join the non-overlapping start of read 1 with read 2.
            r1_insert_start = r1_prefix_length if trim_suffix else r1_overlap_start
            r2_insert_end = len(r2_rc_seq) - r2_prefix_length if trim_suffix else len(r2_rc_seq)
        insert_seq = r1_seq[r1_insert_start: -overlap_size] + r2_rc_seq[: r2_insert_end]
        insert_score = (r1_score[r1_insert_start: -overlap_size]
                        + ''.join(reversed([char for char in r2_score]))[: r2_insert_end])
        r1_precedes_r2 = True
    # MISMATCHES
    else:
        insert_seq = ''
        insert_score = ''
    return insert_seq, insert_score, overlap_size, r1_precedes_r2


def find_longest_read_overlap(seq1, seq2, min_overlap_size=0, seq2_can_precede_seq1=False):
    # Return overlap start index in seq1, overlap start index in seq2, and length of overlap.
    l = len(seq1)
    # For now, reads are the same length.
    assert len(seq2) == l

    # Check full overlap.
    if seq2_can_precede_seq1:
        if seq1 == seq2:
            return 0, 0, l

    shift = 1
    while l - shift >= min_overlap_size:
        # Shift seq1 forward relative to seq2.
        if seq2_can_precede_seq1:
            if seq1[0: -shift] == seq2[shift: l]:
                return 0, shift, l - shift
        # Shift seq1 backward relative to seq2.
        if seq1[shift: l] == seq2[0: -shift]:
            return shift, 0, l - shift
        shift += 1
    return 0, 0, 0


### PYTEST TESTS
def test_merge_read1_read2_fully_overlapping_with_prefixes():
    r1_seq = 'AAAAACCCCCGGGGGTTTTT'
    r2_seq = 'CCCCCGGGGGTTTTTAAAAA' # rc = TTTTTAAAAACCCCCGGGGG
    r1_score = 'EEEEEEEEEEEEEEEJJJJJ'
    r2_score = 'EEEEEEEEEEEEEEEJJJJJ'
    r1_prefix_seq = 'AAA'
    r2_prefix_seq = 'CCC'

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_score,
        r2_score,
        len(r1_prefix_seq),
        len(r2_prefix_seq),
        check_full_overlap=True,
        trim_suffix=True) == ('AACCCCCGG', 'EEEEEEEEE', 15, False)


def test_merge_read1_read2_fully_overlapping_without_prefixes():
    r1_seq = 'AAAAACCCCCGGGGGTTTTT'
    r2_seq = 'CCCCCGGGGGTTTTTAAAAA' # rc = TTTTTAAAAACCCCCGGGGG
    r1_score = 'EEEEEEEEEEEEEEEJJJJJ'
    r2_score = 'EEEEEEEEEEEEEEEJJJJJ'
    r1_prefix_seq = ''
    r2_prefix_seq = ''

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_score,
        r2_score,
        len(r1_prefix_seq),
        len(r2_prefix_seq),
        check_full_overlap=True,
        trim_suffix=True) == ('AAAAACCCCCGGGGG', 'EEEEEEEEEEEEEEE', 15, False)


def test_merge_read1_read2_partially_overlapping_with_prefixes():
    r1_seq = 'AAAAATTTTTGGGGGCCCCC'
    r2_seq = 'TTTTTGGGGGCCCCCAAAAA' # rc = TTTTTGGGGGCCCCCAAAAA
    r1_score = 'JJJJJEEEEEEEEEEEEEEE'
    r2_score = 'JJJJJEEEEEEEEEEEEEEE'
    r1_prefix_seq = 'AAA'
    r2_prefix_seq = 'TTT'

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_score,
        r2_score,
        len(r1_prefix_seq),
        len(r2_prefix_seq),
        check_full_overlap=True,
        trim_suffix=True) == ('AATTTTTGGGGGCCCCCAA', 'JJEEEEEEEEEEEEEEEJJ', 15, True)


def test_merge_read1_read2_partially_overlapping_without_prefixes():
    r1_seq = 'AAAAATTTTTGGGGGCCCCC'
    r2_seq = 'TTTTTGGGGGCCCCCAAAAA' # rc = TTTTTGGGGGCCCCCAAAAA
    r1_score = 'JJJJJEEEEEEEEEEEEEEE'
    r2_score = 'JJJJJEEEEEEEEEEEEEEE'
    r1_prefix_seq = ''
    r2_prefix_seq = ''

    assert merge_read1_read2(
        r1_seq,
        r2_seq,
        r1_score,
        r2_score,
        len(r1_prefix_seq),
        len(r2_prefix_seq),
        check_full_overlap=True,
        trim_suffix=True) == ('AAAAATTTTTGGGGGCCCCCAAAAA', 'JJJJJEEEEEEEEEEEEEEEJJJJJ', 15, True)