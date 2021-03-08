#! /usr/bin/env python3

"""
Read methylation file and print out only rows for requested list of samples.
Handles any file with sample identifier in first column
"""

import sys
import gzip
import matrix_utils

if len(sys.argv) < 3:
    print('\nusage: extract_methyl_samples.py samples.txt methyl.tsv.gz\n')
    sys.exit(1)
    
sample_file = sys.argv[1]
methyl_file = sys.argv[2]

samples = set()

with open(sample_file, 'r') as f:
    for line in f:
        samples.add(line.rstrip())

f = matrix_utils.open_file(methyl_file)
for line in f:
    fields = line.rstrip().split('\t')
    if line.startswith('case') or fields[0] in samples:
        print(line, end='')
f.close()            