# -*- coding: utf-8 -*-

# Copyright (C) 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.


import os
import textwrap

from IlluminaUtils.utils.helperfunctions import remove_spaces

E = os.path.exists
J = os.path.join


class RunConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = remove_spaces(e)
        return
    def __str__(self):
        error_type = 'Config File Error'

        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        error_lines = ['\033[0;30m\033[46m%s%s\033[0m' % (l, ' ' * (max_len - len(l)))\
                                         for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s: %s' % (error_type, error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(error_type) + 2), error_line))

        return '\n'.join(error_message)


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

class RunConfiguration:
    def __init__(self, config):
        self.sanity_check(config)

        self.project_name     = config.get('general', 'project_name').strip()
        self.researcher_email = config.get('general', 'researcher_email').strip()
        self.input_directory  = config.get('general', 'input_directory').strip()
        self.output_directory = config.get('general', 'output_directory').strip()

        self.pair_1 = [os.path.join(self.input_directory, f.strip()) for f in config.get('files', 'pair_1').split(',')]

        if config.has_option('files', 'pair_2'):
            self.pair_2 = [os.path.join(self.input_directory, f.strip()) for f in config.get('files', 'pair_2').split(',')]

        self.trim_to = int(config.get('execute', 'trim_to')) if config.has_option('execute', 'trim_to') else None
        self.min_base_q = int(config.get('execute', 'min_base_q')) if config.has_option('execute', 'min_base_q') else None
        self.ignore_bases = [int(l) for l in config.get('execute', 'ignore_bases').split(',')] if config.has_option('execute', 'ignore_bases') else None
        self.eliminate_Ns = True if (config.get('execute', 'eliminate_Ns') if config.has_option('execute', 'eliminate_Ns') else None) == 'True' else None

        self.pair_1_prefix = config.get('prefixes', 'pair_1_prefix') if config.has_option('prefixes', 'pair_1_prefix') else None
        self.pair_2_prefix = config.get('prefixes', 'pair_2_prefix') if config.has_option('prefixes', 'pair_2_prefix') else None



    def sanity_check(self, config):
        config_template = {
            'general': {
                        'project_name'    : {'mandatory': True},
                        'researcher_email': {'mandatory': True},
                        'input_directory': {'test': lambda x: os.path.exists(x), 'mandatory': True},
                        'output_directory': {'test': lambda x: os.path.exists(x), 'mandatory': True},
            },

            'files': {
                        'pair_1': {'test': lambda x: False not in [E(J(config.get('general', 'input_directory').strip(), t.strip())) for t in x.split(',')], 'mandatory': True,
                                    'required': 'File to be found at %s' % (config.get('general', 'input_directory'))},
                        'pair_2': {'test': lambda x: False not in [E(J(config.get('general', 'input_directory').strip(), t.strip())) for t in x.split(',')],
                                    'required': 'File to be found at %s' % (config.get('general', 'input_directory'))},
            },

            'prefixes': {
                        'pair_1_prefix': {'test': lambda x: len(x) > 0},
                        'pair_2_prefix': {'test': lambda x: len(x) > 0},
            },


            'execute': {
                        'trim_to': {'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 100,
                                    'required': 'Integer value between 1 and 100'},
                        'min_base_q': {'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 40,
                                       'required': 'Integer value between 1 and 40'},
                        'ignore_bases': {'test': lambda x: False not in [RepresentsInt(t) and int(t) > 0 and int(t) <= 100 for t in x.split(',')],
                                         'required': 'Comma separated integers between 1 and 100'},
                        'eliminate_ns': {'test': lambda x: x in ['True', 'False'],
                                         'required': 'True or False'},
            }
        }

        if not len(config.sections()):
            raise RunConfigError('RunConfiguration class is upset, because the config object is empty.\
                                   Probably your config file is not where you think it is. Please check\
                                   your paths.')

        for section in config.sections():
            if section not in config_template:
                raise RunConfigError('Unknown section: "%s"' % (section))
            for option, value in config.items(section):
                if option not in list(config_template[section].keys()):
                    raise RunConfigError('Unknown option under "%s" section: "%s"' % (section, option))
                if 'test' in config_template[section][option] and not config_template[section][option]['test'](value):
                    if 'required' in config_template[section][option]:
                        r = config_template[section][option]['required']
                        raise RunConfigError('Unexpected value for "%s" section "%s": %s \n    Expected: %s' % (option, section, value, r))
                    else:
                        raise RunConfigError('Unexpected value for "%s" section "%s": %s' % (option, section, value))

        for section in config_template:
            for option in config_template[section]:
                if 'mandatory' in config_template[section][option] and not config.has_option(section, option):
                    raise RunConfigError('Missing mandatory option for section "%s": %s' % (section, option))


        if config.has_option('files', 'pair_2'):
            p1 = config.get('files', 'pair_1')
            p2 = config.get('files', 'pair_2')
            if len(p1.split(',')) != len(p2.split(',')):
                raise RunConfigError('For paired end setup, number of files has to match in pair_1 and pair_2.')
