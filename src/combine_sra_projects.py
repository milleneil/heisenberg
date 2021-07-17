#! /usr/bin/env python3

""" 
combine list of methyl 'wide' files into one and ensure that all probe vals are
in the same column order for each sample
"""

import sys
import gzip
import argparse
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='safely combine multiple files',
                                     prog="heisenberg combine")
    parser.add_argument('-i', '--input', type=str, 
                        help='source file of methylation values', 
                        required=True, action='append')

    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')

    parser.add_argument('-p', '--probes', type=str,
                        help='subset of probes to ensure are in output file ' + 
                        '[default=inferred from 1st line of 1st file]')
    
    parser.add_argument('-x', '--probe_start_idx', type=int, default=4,
                        help='start index of methyl probes [default=4]')
    
    parser.add_argument('-m', '--missing_val', type=float, default=-1,
                        help='placeholder value to use for missing probe vals [default=-1]')
    
        
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    probes = []
    # read probes from file if supplied, else read from input file
    if args.probes:
        with open(args.probes, 'r') as f:
            for line in f:
                probes.append(line.rstrip())
    else:
        f = matrix_utils.open_file(args.input[0])
        header = f.readline().rstrip().split('\t')
        probes = header[args.probe_start_idx:]
        f.close()
      
    # sort to put probes in order roughly
    probes = sorted(probes)       
    print(f'{len(probes)} probes loaded', file=sys.stderr)
 
    out_file = sys.stdout
    if args.output:
        out_file = matrix_utils.open_output_file(args.output)   

    header_cols = ['case','sample','biospecimen','tissue']
    header_cols.extend(probes)
    print('\t'.join(header_cols), file=out_file)
    
    for file in args.input: 
        print(f"reading {file}", file=sys.stderr)
        f = matrix_utils.open_file(file)
        # read header and record col index for each probe name
        header = f.readline().rstrip().split('\t')
        header_map = {}
        for i in range(0,len(header)):
            header_map[header[i]] = i
            
        col_idxs = []
        for p in probes:
            col_idxs.append(header_map[p])
                
        for line in f:
            fields = line.rstrip().split('\t')
            vals = fields[0:4]
            for idx in col_idxs:
                vals.append(fields[idx])
            print('\t'.join(vals), file=out_file)
        f.close()
    out_file.close()
            
if __name__ == "__main__":
    main()    