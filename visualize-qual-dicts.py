# -*- coding: utf-8 -*-

import os
import sys
import cPickle

sys.path.append('/bioware/seqinfo/illumina-utils/')
import fastqlib as u

def main(qual_dict, figure_dest, title = None, split_tiles = False):
    plot_dict = u.compute_plot_dict_from_tiles_dict(qual_dict)
    u.visualize_qual_stats_dict(plot_dict, figure_dest, title, split_tiles)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Visualize Qual Dicts')
    parser.add_argument('qual_dict', metavar = 'QUAL_DICT',
                                        help = 'cPickle dictionary that contains quality score info')
    parser.add_argument('-d', '--dest', metavar = 'DEST_FILE', default = None,
                                        help = 'Figure destination')
    parser.add_argument('--title', metavar = 'TITLE', default = None,
                                        help = 'Title to appear at the top of the figure')
    parser.add_argument('--split-tiles', action = 'store_true', default = False,
                                        help = 'When set, quality curves will be shown separately for each tile')

    args = parser.parse_args()
    
    qual_dict = cPickle.load(open(args.qual_dict))
    dest = args.dest if args.dest else args.qual_dict
    title = args.title if args.title else os.path.basename(args.qual_dict)
    split_tiles = args.split_tiles

    main(qual_dict, dest, title, split_tiles)

