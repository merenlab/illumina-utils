# -*- coding: utf-8 -*-

import os
import sys
import ConfigParser

E = os.path.exists
J = os.path.join

import fastqlib as u

class ConfigError(Exception):
    pass

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

class RunConfiguration:
    def __init__(self, config):
        self.sanity_check(config)

        self.researcher_email = config.get('general', 'researcher_email').strip()
        self.input_directory  = config.get('general', 'input_directory').strip()
        self.output_directory = config.get('general', 'output_directory').strip()

        self.lane_1 = [os.path.join(self.input_directory, f.strip()) for f in config.get('files', 'lane_1').split(',')]

        if config.has_option('files', 'lane_2'):
            self.lane_2 = [os.path.join(self.input_directory, f.strip()) for f in config.get('files', 'lane_2').split(',')]

        self.trim_to = int(config.get('execute', 'trim_to')) if config.has_option('execute', 'trim_to') else None
        self.min_base_q = int(config.get('execute', 'min_base_q')) if config.has_option('execute', 'min_base_q') else None
        self.ignore_bases = [int(l) for l in config.get('execute', 'ignore_bases').split(',')] if config.has_option('execute', 'ignore_bases') else None
        self.eliminate_Ns = True if (config.get('execute', 'eliminate_Ns') if config.has_option('execute', 'eliminate_Ns') else None) == 'True' else None


        
    def sanity_check(self, config):
        config_template = {
            'general': {
                        'researcher_email': {'mandatory': True},
                        'input_directory': {'test': lambda x: os.path.exists(x), 'mandatory': True},  
                        'output_directory': {'test': lambda x: os.path.exists(x), 'mandatory': True},
            },

            'files': {
                        'lane_1': {'test': lambda x: False not in [E(J(config.get('general', 'input_directory').strip(), t.strip())) for t in x.split(',')], 'mandatory': True},
                        'lane_2': {'test': lambda x: False not in [E(J(config.get('general', 'input_directory').strip(), t.strip())) for t in x.split(',')]},
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

        for section in config.sections():
            if section not in config_template:
                raise ConfigError, 'Unknown section: "%s"' % (section)
            for option, value in config.items(section):
                if option not in config_template[section].keys():
                    raise ConfigError, 'Unknown option under "%s" section: "%s"' % (section, option)
                if config_template[section][option].has_key('test') and not config_template[section][option]['test'](value):
                    if config_template[section][option].has_key('required'):
                        r = config_template[section][option]['required']
                        raise ConfigError, 'Unexpected value for "%s" section "%s": %s \n    Expected: %s' % (option, section, value, r)
                    else:
                        raise ConfigError, 'Unexpected value for "%s" section "%s": %s' % (option, section, value)

        for section in config_template:
            for option in config_template[section]:
                if config_template[section][option].has_key('mandatory') and not config.has_option(section, option):
                    raise ConfigError, 'Missing mandatory option for section "%s": %s' % (section, option)


        if config.has_option('files', 'lane_2'):
            l1 = config.get('files', 'lane_1')
            l2 = config.get('files', 'lane_2')
            if len(l1.split(',')) != len(l2.split(',')):
                raise ConfigError, 'For paired end setup, number of files has to match in lane_1 and lane_2.'
                    


def main(config):
    #
    # FIXME
    #
    # At this point, this is not the ideal main function at all. These are what I need:
    #   - define an ordered list of filter functions for every sequence to go through,
    #   - whatever filter is defined in config.execute, create a reduced set of filters,
    #   - call them like a boss.
    #
    # but right now I am very lazy to start that part, so I will only implement eliminate_ns
    # filter right below. but it needs to be carried in a module very soon.
    #

    errorstxt = open(os.path.join(config.output_directory, 'error_counts.txt'), 'w')

    for index in range(0, len(config.lane_1)):
        try:
            input_1 = u.FastQSource(config.lane_1[index], compressed = True)
            input_2 = u.FastQSource(config.lane_2[index], compressed = True)

            output_1 = u.FastQOutput(os.path.join(config.output_directory, os.path.basename(config.lane_1[index][:-3])))
            output_2 = u.FastQOutput(os.path.join(config.output_directory, os.path.basename(config.lane_2[index][:-3])))
        except u.FastQLibError, e:
            print "FastQLib is not happy.\n\n\t", e, "\n"
            sys.exit()

        total_number_of_reads = 0
        ns_in_pair_2 = 0
        ns_in_pair_1 = 0
        
        while input_1.next() and input_2.next():
            if input_1.p_available:
                input_1.print_percentage()
            total_number_of_reads += 1
            if 'N' in input_1.entry.sequence:
                ns_in_pair_1 += 1
            elif 'N' in input_2.entry.sequence:
                ns_in_pair_2 += 1
            else:
                output_1.store_entry(input_1.entry)
                output_2.store_entry(input_2.entry)

        print
        number_of_pairs_with_ns = ns_in_pair_1 + ns_in_pair_2
        errorstxt.write('file name              : %s\ntotal number of reads  : %d\nnumber of pairs with ns: %d (%%%.2f)\nns in pair 1           : %d (%%%.2f)\nns in pair 2           : %d (%%%.2f)\n\n' % (os.path.basename(config.lane_1[index]), total_number_of_reads, number_of_pairs_with_ns, number_of_pairs_with_ns * 100.0 / total_number_of_reads, ns_in_pair_1, ns_in_pair_1 * 100.0 / number_of_pairs_with_ns, ns_in_pair_2, ns_in_pair_2 * 100.0 / number_of_pairs_with_ns))
    errorstxt.close()
    print


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Microbial Pseudo Community Sample Generation')
    parser.add_argument('user_config', metavar = 'CONFIG_FILE',
                                        help = 'User configuration to run')

    args = parser.parse_args()
    user_config = ConfigParser.ConfigParser()
    user_config.read(args.user_config)

    try: 
        config = RunConfiguration(user_config)
    except ConfigError, e:
        print "There is something wrong with the config file. This is what we know: \n\n", e
        print
        sys.exit()

    sys.exit(main(config))
