#! /usr/bin/env python3

import statistics
import argparse
import matrix_utils
import sys

"""
Read input .wide.tsv probes file with probes as columns and calculate
descriptive statistics for each probe
"""


def parse_args():
    parser = argparse.ArgumentParser(description='gather descriptive statistics' +
                                                ' for all probes in input',
                                     prog="heisenberg stats")
    parser.add_argument('-i', '--input', type=str, 
                        help='source file of methylation values', required=True)

    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')

    parser.add_argument('-p', '--probe_start_idx', type=int, default=4,
                        help='column index of first methyl probe in input [default=4]')
    
    parser.add_argument('-m', '--missing_val', type=float, default=-1,
                        help='placeholder value to use for missing probe vals [default=-1]')
    
    parser.add_argument('-x', '--max_probes', type=float, default=25000,
                        help='number of probes to process at once [default=25k]')

        
    args = parser.parse_args()
    return args

def read_header(in_file, probe_start_idx):

    print(f'reading header from {in_file}...', file=sys.stderr)
    probe_idxs = {}
    f = matrix_utils.open_file(in_file)
    for line in f:
        fields = line.rstrip().split('\t')
        for i in range(probe_start_idx, len(fields)):
            probe_idxs[fields[i]] = i
        break
    f.close()
    return probe_idxs    

def split_probes(probe_idxs, chunk_size):
    probe_chunks = []
    current_chunk = {}
    for probe, _ in probe_idxs.items():
        if len(current_chunk) == chunk_size:
            probe_chunks.append(current_chunk)
            current_chunk = {}
        current_chunk[probe] = []
    probe_chunks.append(current_chunk)
    return probe_chunks
            

def main():
    """
	Read input file and save values from each row for each probe. After reading
	file, calculate statistics and print.
	
	This script reads all rows for target probe set and temporarily stores
	values in memory - full 450k requires > 16GB RAM - by default, file 
	is processed repeatedly in chunks of 25k.
    """

    args  = parse_args()

    probe_idxs = read_header(args.input, args.probe_start_idx)
    print(f'{len(probe_idxs)} probe indexes loaded', file=sys.stderr)

    output = sys.stdout
    if args.output:
        output = matrix_utils.open_output_file(args.output)

    
    # read probes in multiple passes over file to handle a memory friendly
    # amount of data each time
    probe_chunks = split_probes(probe_idxs, args.max_probes)
    print(f'{len(probe_chunks)} probe splits created', file=sys.stderr)
    chunk_counter = 1
    for probe_vals in probe_chunks :
        print(f'reading probe chunk {chunk_counter} from {args.input}...', 
              file=sys.stderr)
        
        chunk_counter += 1
        
        f = matrix_utils.open_file(args.input)
        # read past header
        f.readline()
        
        for line in f:
            fields = line.rstrip().split('\t')
            
            for probe in probe_vals.keys():
                probe_idx = probe_idxs[probe]
                probe_val = None
                try:
                    probe_val = float(fields[probe_idx])
                    probe_vals[probe].append(probe_val)
                except ValueError:
                    if args.missing_val:
                        probe_vals[probe].append(args.missing_val)
        f.close()
                            
        print(f'calculating stats and writing to output...', file=sys.stderr)
    
        print('\t'.join(['probe','min','max','mean','stdev']), file=output)
        for probe, vals in probe_vals.items():
            minimum = min(vals)
            maximum = max(vals)
            mean = statistics.mean(vals)
            stdev = statistics.stdev(vals)
            printvals = [probe,str(minimum),str(maximum),str(mean),str(stdev)]
            print('\t'.join(printvals), file=output)
        # free up mem
        probe_vals.clear()
    output.close()
    print('completed', file=sys.stderr)

if __name__ == "__main__":
    main()    
