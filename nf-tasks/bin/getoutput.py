#!/usr/bin/env python3
#for: humann.nf

import sys, os
import random

scriptSelfDirectory = os.path.dirname(os.path.realpath(__file__))
dataZip = os.path.join(scriptSelfDirectory, "DATA/samples.zip")
templates = ('Sample1', 'Sample2', 'Sample3')
# Get argument from command line
if len(sys.argv) < 2:
    exit()
sampleName = sys.argv[1]

def extractSample(name, zipfile, prefix='Sample'):
    # name: one of 'Sample1', 'Sample2', 'Sample3' (see "templates")
    # zipfile: path to the zip file
    # prefix: prefix of the sample file name

    files = (f'{name}_genefamilies.tsv', f'{name}_metaphlan_bugs_list.tsv', f'{name}_pathabundance.tsv', f'{name}_pathcoverage.tsv')
    for file in files:
        # output is file replacing name with prefix
        output = file.replace(name, prefix)
        os.system(f'unzip -p {zipfile} {file} > {output}')
    
        # Edit the first line of {output} replacing with regex /(\S+)_Abundance-RPKs/ with /{name}_Abundance-RPKs/
        keystrings = ('_Abundance-RPKs', '_Coverage', '_Abundance')
        with open(f'{output}', 'r') as f:
            lines = f.readlines()
            for index, line in enumerate(lines):
                parts = line.split('\t')

                # if any of the keystrings is in parts[1]

                if len(parts) == 2 and any(key in parts[1] for key in keystrings):
                    (sample, suffix) = parts[1].split('_')
                    parts[1] = f'{prefix}_{suffix}'
                
                    lines[index] = '\t'.join(parts)
                    #print(f"Found keystring {parts[1]} at line {index}: {lines[index]}", file=sys.stderr)
                    with open(f'{output}', 'w') as f:
                        f.writelines(lines)
                        break



if sampleName in templates:
    template = sampleName
else:
    # pick a random sample
    template = random.choice(templates)

print(f"Extracting '{template}' as '{sampleName}'" , file=sys.stderr)
extractSample(template, dataZip, prefix=sampleName)