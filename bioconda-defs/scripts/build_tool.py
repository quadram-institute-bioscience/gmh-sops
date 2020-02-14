#!/usr/bin/env python3

import sys
import argparse
import pprint
import yaml
import subprocess
import os

def eprint(*args, **kwargs):
    """print to STDERR"""
    print(*args, file=sys.stderr, **kwargs)

def vprint(*args, **kwargs):
    if not opt.verbose:
        return 0
    eprint(*args, **kwargs)

def save_to_file(string, output):
	eprint("\tSave to {}".format(output))
	if os.path.isfile(output):
		eprint(" *** WARNING: Skipping {}, file was found".format(output))
		return False
	try:
		file = open(output, 'w')
		file.write(string.decode("utf-8"))
		file.close()
	except Exception as e:
		eprint("ERROR:\nUnable to write to '{}': {}".format(output, e))
		exit(1)

def make_definition(package, entrypoint):
	eprint("\tMake def for {}".format(package))
	# Clean inputs
	command = ['python3', defscript, '-c', 'bioconda' , '-p', package, '-x', entrypoint.rstrip()]
	try:
		cmd = ['python3', defscript, '-c', 'bioconda' , '-p', package, '-x', entrypoint.rstrip()]
		definition = subprocess.check_output(cmd)
		return definition
	except Exception as e:
		eprint("ERROR:\nUnable to execute external script: {}.\n{}".format(e, cmd))
		exit(1)


defscript = os.path.dirname(os.path.realpath(__file__)) + '/conda2def.py'

if len(sys.argv)<2:
	eprint("Missing argument: list file in PACKAGE=VERSION{tab}ENTRYPOINT format")
else:
	with open(sys.argv[1]) as f:
		lines = f.readlines()
		for line in lines:
			[package, entrypoint] = line.split("\t")
			[package_name, package_version] = package.split('=')
			eprint("[{}] version={};cmd={};".format(package_name, package_version, entrypoint.rstrip() ) )
			outfile = os.path.dirname(os.path.realpath(sys.argv[1])) + '/' + package_name + '-' + package_version + '.def'
			definition = make_definition(package, entrypoint)
			save_to_file(definition, outfile)
