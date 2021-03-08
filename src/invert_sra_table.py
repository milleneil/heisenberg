#! /usr/bin/env python3

"""
read file of methylation beta values with probes as rows and samples as
columns and turn into file with samples and rows and probes as columns

print extra columns for sample/biospecimen and tissue type to match
format for TCGA wide methylation files -- tissue type for all samples
is peripheral blood
"""

import sys
import gzip
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='invert table with probes ' +
                                     'as rows to probes as columns',
                                     prog="heisenberg invert")
    parser.add_argument('-i', '--input', type=str, 
                        help='input file with probe values as rows', 
                        required=True)

    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')

    parser.add_argument('-t', '--tissue', type=str, 
                        help='value to write in tissue column [default=normal ' +
                        'peripheral blood]', default='normal peripheral blood')

        
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    biospecimen_probes = {}
    biospecimen_indexes = None
    probes = []
    f = matrix_utils.open_file(args.input)
    
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if fields[0].startswith('ID_REF'):
            biospecimen_indexes = fields
            for i in range(1,len(fields)):
                biospecimen_probes[fields[i]] = {}
        else:        
            probe_id = fields[0]
            probes.append(probe_id)
            for i in range(1, len(biospecimen_indexes)):
                if i >= len(fields):
                    probe_val = 'NA'
                else:
                    probe_val = fields[i]
                    if probe_val == '':
                            probe_val = 'NA'
                biospecimen = biospecimen_indexes[i]
                biospecimen_probes[biospecimen][probe_id] = probe_val
    f.close()

    output = sys.stdout
    if args.output:
        output = matrix_utils.open_output_file(args.output)

    headers = ['case', 'sample', 'biospecimen', 'tissue']
    for probe in probes:
        headers.append(probe)
    print("\t".join(headers), file=output)
    for biospecimen, probe_vals in biospecimen_probes.items():
        vals = []
            # same identifiers x3 to match other files
        vals.append(biospecimen)
        vals.append(biospecimen)
        vals.append(biospecimen)
        vals.append(args.tissue)
        for probe in probes:
            probe_val = probe_vals[probe]
            vals.append(probe_val)
        print("\t".join(vals), file=output)
    
    output.close()
    
if __name__ == "__main__":
    main()    