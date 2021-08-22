#! /usr/bin/env python3

# print specific cell value of a matrix file given row/col labels

import sys
import gzip
import argparse
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='print value of matrix cell',
                                     prog="heisenberg cell")
    parser.add_argument('-i', '--input', type=str, 
                        help='source file of methylation values', 
                        required=True)

    parser.add_argument('-c', '--column', type=str, 
                        help='column label to print', required=True)
    parser.add_argument('-r', '--row', type=str, help='row label to print',
                        required=True)
        
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    # assume first line of file is headers and first column is labels
    f = matrix_utils.open_file(args.input)
        
    col_idx = None
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if col_idx == None:
            for i in range(0,len(fields)):
                if fields[i] == args.column:
                    col_idx = i
                    print(f"col idx: {col_idx} for {args.column}", file=sys.stderr)
                    break
            if col_idx == None:
                raise Exception("col header not found: " + args.column)
        else:
            if fields[0] == args.row:
                print(f"\t{args.column}")
                print(f"{args.row}\t{fields[col_idx]}")
                break
            
if __name__ == "__main__":
    main()    