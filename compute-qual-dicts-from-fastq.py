# -*- coding: utf-8 -*-

import os
import sys
import numpy
import cPickle
import ConfigParser

sys.path.append('/bioware/seqinfo/illumina-utils/')

import fastqlib as u
from runconfiguration import RunConfiguration
from runconfiguration import ConfigError


def main(input_1_path, input_2_path, output_prefix):
   
    ########################################################################################################################
    # qual dicts
    ########################################################################################################################
    
    tiles_dict  = {'1': {}, '2': {}}

    input_1 = u.FastQSource(input_1_path)
    input_2 = u.FastQSource(input_2_path)
    
    tiles_dict = u.populate_tiles_qual_dict_from_input(input_1, input_2, tiles_dict)

    sys.stderr.write('[qual dicts] computing plot dict ... ')
    plot_dict = u.compute_plot_dict_from_tiles_dict(tiles_dict)
    sys.stderr.write('done.\n')

    sys.stderr.write('[qual dicts] serializing plot dicts ... ')
    cPickle.dump(plot_dict, open(os.path.join(output_prefix), 'w'))
    sys.stderr.write('done.\n')
 
    input_1.close()
    input_2.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
