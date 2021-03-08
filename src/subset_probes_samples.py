#! /usr/bin/env python3

# extract subset of probe values from full input set

import sys
import argparse
import matrix_utils

def parse_args():
    parser = argparse.ArgumentParser(description='extract methylation values ' +
                                     'for probe and/or sample subset',
                                     prog='heisenberg subset')
    
    parser.add_argument('-i', '--input', type=str, 
                        help='source file of methylation values', required=True)
    
    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')
    
    parser.add_argument('-s', '--samples', type=str, help='samples to extract')

    parser.add_argument('-p', '--probes', type=str, help='probess to extract')
    
    parser.add_argument('-x', '--probe_start_idx', type=int, default=4,
                        help='column index of first methyl probe in input [default=4]')
    
    parser.add_argument('-l', '--lenient', action="store_true",
                        help='tolerate missing probes')
    
    
    args = parser.parse_args()
    
    return args

def load_file(file):
    """ load sample or probe list """
    vals = set()
    with open(args.samples, 'r') as f:
        for line in f:
            # skip header if found
            if line.startswith('Illumina'):
                continue
            # accommodate case where sample file is a simple list of names 
            # as well as when it's a delimited file (but assume sample is 
            # first col)
            fields = line.rstrip('\n').split('\t')
            vals.add(fields[0])
    return vals

def main():
    args = parse_args()
    
    samples = None
    probes = None
    
    if args.samples:
        print(f'reading samples from: {args.probes}', file=sys.stderr)
        samples = load_file(args.samples)
        print(f'{len(samples)} samples loaded', file=sys.stderr)
         
    if args.probes:
        print(f'reading probes from: {args.probes}', file=sys.stderr)
        probes = load_file(args.probes)
        print(f'{len(probes)} probes loaded', file=sys.stderr)
    
    output = sys.stdout
    if args.output:
        print(f'writing to output file: {args.output}', file=sys.stderr)
        output = matrix_utils.open_output_file(args.output)
    
    probe_idxs = []
    found_probes = set()
    header = True
    
    print(f'reading from input file: {args.input}', file=sys.stderr)
    f = matrix_utils.open_file(args.input)
    
    counter = 0
    # default to sample index of 1
    sample_idx = 1
    for line in f:
        counter += 1
        if counter % 10000 == 0:
            print(f'{counter} lines read from input', file=sys.stderr)
        
        fields = line.rstrip().split('\t')
        if header:
            # read column headers and identify column indexes for probe subset
            for i in range(0, len(fields)):

                if fields[i] == 'sample':
                    sample_idx = i
                
                # keep metadata columns of sample info no matter what as 
                # well as cols in our probe subset
                if probes == None or i < args.probe_start_idx or fields[i] in probes:
                    probe_idxs.append(i)
                    found_probes.add(fields[i])
                
            # make sure we got all probes
            if probes != None:
                for probe in probes:
                    if not probe in found_probes:
                        msg = "probe not found in header: " + probe
                        if args.lenient:
                            print(msg, file=sys.stderr)
                        else:
                            raise Exception(msg)
            
            header = False
        
        vals = []
        # filter by samples if requested
        if samples == None or fields[sample_idx] in samples:
            for idx in probe_idxs:
                vals.append(fields[idx])
            print('\t'.join(vals), file=output)
    
    output.close()
    f.close()    
    
    print(f'{counter} lines written...completed', file=sys.stderr)

if __name__ == "__main__":
    main()    
