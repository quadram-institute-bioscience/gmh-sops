#!/usr/bin/env python3
"""
Create a singularity definition file from a list of conda packages
or from an ENVIRONMENT.YAML file.
"""

import sys
import argparse
import pprint
import yaml

def eprint(*args, **kwargs):
    """print to STDERR"""
    print(*args, file=sys.stderr, **kwargs)

def vprint(*args, **kwargs):
    if not opt.verbose:
        return 0 
    eprint(*args, **kwargs)

def makeDefFromYaml():
    print("""
Bootstrap: docker

From: continuumio/miniconda3

%environment
 PATH=/opt/conda/envs/{1}/bin:$PATH

%post
 echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
 echo "source activate {1}" > ~/.bashrc
 /opt/conda/bin/conda create -n {1} {2} -y {3}

%runscript
 exec {4} "$@"
""".format(
    opt.env_file,
    env_name,
    channels_string,
    packages_string,
    opt.cmd
))
    

def makeDefFromList():
    pass

opt_parser = argparse.ArgumentParser(description='Create a Singularity definition file from a list of Conda packages')

opt_parser.add_argument('-e', '--env-file',
                        help='Conda environment file in YAML format')

opt_parser.add_argument('-c', '--channels',
                        help='Conda channel (when not using --env-file)',
                        action='append')

opt_parser.add_argument('-p', '--packages',
                        help='Conda package(s) (when not using --env-file)',
                        action='append')

opt_parser.add_argument('-b', '--ignore-build',
                        help='Ignore package build number (use when YAML file was generated in OSX)',
                        action='store_true')

opt_parser.add_argument('-x', '--cmd',
                        help='Default command',
                        default='')


opt_parser.add_argument('-r', '--ignore-version',
                        help='Ignore package version (not recommended, implies -b)',
                        action='store_true')

opt_parser.add_argument('-v', '--verbose',
                        help='Increase output verbosity',
                        action='store_true')


opt = opt_parser.parse_args()


if __name__ == '__main__':
    
    if opt.env_file != None:
        vprint('- Environment mode: {}'.format(opt.env_file))
        try:
            with open(opt.env_file, 'r') as stream:
                try:
                    data = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    eprint('\nFATAL YAML ERROR:\ntrying to parse <{}>:\n{}'.format(opt.env_file,exc))
        except Exception as exc:
                    eprint('\nFATAL ERROR:\nTrying to read {}:\n{}'.format(opt.env_file,exc))

        env_name = 'default'
        if 'name' in data:
            env_name = data['name']
        channels_string = ''
        packages_string = ''
        if data['channels'] != None:
            for c in data['channels']:
                channels_string += ' -c {} '.format(c)

        for p in data['dependencies']:
            package_list = p.split('=')
            if len(package_list) < 2:
                eprint("ERROR: Package {} is not in the 'name'='version'='build' format".format(p))
                exit(1)
            if (opt.ignore_version):
                packages_string += '{} '.format(package_list[0])
            elif (opt.ignore_build):
                packages_string += '{}={} '.format(package_list[0], package_list[1])
            else:
                packages_string += ' {} '.format(p)

        makeDefFromYaml()
    elif opt.packages != None:
        vprint('- Package mode:')

        env_name = 'container'
        channels_string = ''
        packages_string = ''
        if opt.channels != None:
            for channel in opt.channels:
                vprint(' - Channel: {}'.format(channel))
                channels_string += ' -c {} '.format(channel)

        for package in opt.packages:
            package_list = package.split('=')
            if (opt.ignore_version):
                packages_string += '{} '.format(package_list[0])
                vprint(' - Package (no version): {} (was: {})'.format(package_list[0], package))
            elif (opt.ignore_build) and (len(package_list) > 1):
                packages_string += '{}={} '.format(package_list[0], package_list[1])
                vprint(' - Package (no build): {}={} (was: {})'.format(package_list[0], package_list[1], package))
            else:
                packages_string += ' {} '.format(package)
                vprint(' - Package: {}'.format(package))

        makeDefFromYaml()
    else:
        eprint("Missing parameters: Specify at least a package with -p PACKAGE (or env file)\nUse -h or --help to print usage.")
        exit()

