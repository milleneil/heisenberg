#! /usr/bin/env python3

"""
Read SRA series_matrix.txt.gz file and pull out sample methylation
vals for each probe -- print matrix file with sample ids as columns and
methylation probes as rows
"""

import sys
import gzip
import argparse
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='extract probe values from SRA file',
                                     prog="heisenberg extract_sra_probe")
    parser.add_argument('-i', '--input', type=str, 
                        help='series_matrix.txt file of methylation values', 
                        required=True)

    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')
        
    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    matrix_started = False
    
    output = sys.stdout
    if args.output:
        output = matrix_utils.open_output_file(args.output)
    
    f = matrix_utils.open(args.input)
    for line in f:
        # skip all lines until matrix table starts, then print remainder
        if line.startswith('!series_matrix_table_begin'):
            matrix_started = True
        elif line.startswith('!series_matrix_table_end'):
            matrix_started = False
        elif matrix_started:
            # remove all double quotes around vals then print line
            print(line.rstrip().replace('"',''), file=output)
    f.close()
    output.close()

if __name__ == "__main__":
    main()    
