#! /usr/bin/env python3

"""
Read SRA series_matrix.txt.gz file and pull out sample identifier
and requested Sample_characteristics_ch secondary metadata field

eg. 'gender: ', 'sample group: '
"""

import sys
import gzip
import argparse
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='extract SRA metadata for ' +
                                     'samples',
                                     prog='heisenberg extract_sra_meta')
    parser.add_argument('-f', '--file', type=str, 
                        help='series_matrix.txt.gz file', required=True)

    parser.add_argument('-m', '--meta', type=str, 
                        help='sample metadata field to print', required=True, 
                        action='append')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')
    parser.add_argument('-p', '--prefix', type=str, 
                        help='start each row with this')
    parser.add_argument('-t', '--transform', type=str, action='append',
                        help='transform Y to X in output [Y=X]')
        

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    transforms = {}
    if args.transform:
        for val in args.transform:
            pieces = val.split('=')
            transforms[pieces[0]] = pieces[1]
    samples = None
    vals_map = {}
    
    out_file = sys.stdout
    if args.output:
        out_file = open(args.output, 'w')
    
    f = matrix_utils.open(args.file)
    for line in f:
            if line.startswith('!Sample_geo_accession'):
                fields = line.rstrip().replace('"','').split('\t')
                samples = fields[1:]
            else:
                for meta in args.meta:
                    meta_key = '"' + meta + ':'
                    if line.startswith('!' + meta):
                        fields = line.rstrip().replace('"','').split('\t')
                        vals_map[meta] = fields[1:]
                    elif line.startswith('!Sample_characteristics_ch') and meta_key in line:
                        fields = line.rstrip().replace('"','').split('\t')
                        vals_map[meta] = [i.replace(meta + ': ','') for i in fields[1:]]
    f.close()
                            
    
    for i in range(0, len(samples)):
        sample = samples[i]
        vals = [sample]
        if args.prefix:
            vals.insert(0, args.prefix)
        for meta, meta_vals in vals_map.items():
            # apply data transformation if needed
            display_val = meta_vals[i]
            if meta_vals[i] in transforms:
                display_val = transforms[meta_vals[i]]
            
            vals.append(display_val)
        print('\t'.join(vals), file=out_file)
    out_file.close()
        

if __name__ == "__main__":
    main() 
