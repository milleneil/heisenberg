#! /usr/bin/env python3

"""
Make simulated samples by combining probe vals from k chosen samples (i.e.
n choose k where n=total number of input samples).

Choose samples, then create new sample by taking mean of each probe 
value and printing that

Basic formula:

S1 = sample 1 methylation val
S2 = sample 2 methylation val
S3 = simulated sample 3 val
e = random error to add in for noise

S3 = Adjusted methyl value = mean(S1,S2) + e

"""


import sys
import argparse
import matrix_utils
import methyl_sample_utils as m_utils
import simulation_noise
import itertools
import statistics

def parse_args():
    parser = argparse.ArgumentParser(description='create simulated methylation ' +
                                     'probe values',
                                     prog='heisenberg simulate')
    parser.add_argument('-i', '--input', type=str, 
                        help='source file of methylation values', required=True)

    parser.add_argument('-k', '--choose', type=int, 
                        help='number of individuals to choose for simulated sample',
                        default=2)

    parser.add_argument('-x', '--max_individuals', type=int,
                        help='max number of individuals to create in each ' +
                        'age group [default=max possible]')

    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')

    parser.add_argument('-b', '--tissue', type=str, 
                        help='tissue display val for simulated samples ' +
                        '[default=normal peripheral blood]',
                        default='normal peripheral blood')
    
    parser.add_argument('-s', '--stage', type=str, 
                        help='stage display val for simulated samples ' +
                        '[default=none',
                        default='none')
    
    parser.add_argument('-r', '--with_replacement', action='store_true',
                        help='choose samples with replacement [default=false]')
    
    parser.add_argument('-v', '--structural_variants', type=str,
                        help='structural variant probe overlap file')

    parser.add_argument('-c', '--confounding_snps', type=str,
                        help='list of SNPs that may confound methyl values')
    
    parser.add_argument('-a', '--probe_stats', type=str,
                        help='descriptive stats per probe for adding random noise')
    
    parser.add_argument('-l', '--required_probes', type=str,
                        help='full list of probes to be in output regardless of presence in input [pad with -1]')
    
    parser.add_argument('-p', '--probe_start_idx', type=int, default=4,
                        help='column index of first methyl probe in input [default=4]')
    
    parser.add_argument('-y', '--required_only', action='store_true',
                        help='print only required probes (skip all others)')
    
    parser.add_argument('-cd', '--max_distance', type=int, default=2,
                        help='exclude confounding snps greater than this bp away [default=2]')

    parser.add_argument('-cf', '--max_snp_maf', type=float, default=0.2,
                        help='exclude confounding snps with MAF greater than ' +
                            'this [default=.2]')

    parser.add_argument('-sf', '--max_sv_maf', type=float, default=0.3,
                        help='exclude structural variants with MAF greater ' +
                        'than this [default=.3]')
    
    
    args = parser.parse_args()
    return args


def make_combinations(samples, k, with_replacement, out_file, max_samples, tissue, stage, probes):
    """
    Make k unique combinations from samples up to max and print to out_file.
    Return number of combinations created. If max_samples is None, all possible
    combinations will be created
    """
    
    counter = 0
    iterator = None
    if with_replacement:
        iterator = itertools.combinations_with_replacement(samples.values(), k)
    else:
        iterator = itertools.combinations(samples.values(), k)
        
    for combination in iterator:
        if max_samples != None and counter >= max_samples:
            break
        else:
            combined = combine(combination, tissue, stage, probes)
            print(combined, file=out_file)
            counter += 1
    return counter

def aggregate_probe_vals(samples, probes):
    """
    create probe vals that represent mean of each probe across samples
    """
    probe_vals = {}
    
    for sample in samples:
        # assume all samples have the same number of probes and use one at random
        # as model
        
        for probe, val in sample.probe_vals.items():
            if probe not in probe_vals:
                probe_vals[probe] = 0      
            probe_vals[probe] += val
  
    num_samples = len(samples)
    
    adjusted = {}
    
    # for each probe label
        # randomly simulate sv
        # randomly simulate confounding snp
        # simulate random noise based off of stdev
        # adjust aggregate probe val accordingly
    for probe_label, total_val in probe_vals.items():
        avg = total_val / num_samples
        if probe_label in probes:
            probe = probes[probe_label]
            
            # negative values signal missing probe data in one or more
            # samples used to create simulated sample -- skip adjustment
            # in this case and print missing val (-1) for this probe
            if avg > 0:
                avg = simulation_noise.random_adjust(avg, probe)
                adjusted[probe_label] = avg
            else:
                #FIXME: parameterize missing val
                adjusted[probe_label] = -1
        else:
            adjusted[probe_label] = avg
  
    return adjusted

def combine(samples, tissue, stage, probes):
    """ make combined sample with methyl vals in samples"""
    combined = m_utils.MethylSample()        
    combined.case = 'sim'
    combined.sample = '_'.join(set([s.case for s in samples]))
    combined.biospecimen = 'mean-methyl-sim'
    combined.tissue = tissue
    combined.stage = stage
    
    print(f'creating simulated sample {combined.sample}', file=sys.stderr)
    combined.probe_vals = aggregate_probe_vals(samples, probes)
 
    # samples should have same gender and age group, set synthetic sample
    # to be gender and mean of ages
    for sample in samples:
        if hasattr(sample, 'gender'):
            combined.gender = sample.gender
            combined.age = statistics.mean([s.age for s in samples])
            break
    return combined


def main():
    args = parse_args()
    
    # load list of required probes
    required = None
    if args.required_probes:
        print(f'loading required probes from {args.required_probes}', file=sys.stderr)
        required = matrix_utils.load_probe_list(args.required_probes)
        print(f'{len(required)} required probes loaded', file=sys.stderr)
    
    # input files are rows = samples, cols = probe vals
    print(f'reading samples from {args.input}', file=sys.stderr)
    samples = m_utils.load_file_as_samples(args.input, 
                                           required=required, 
                                           start_idx=args.probe_start_idx,
                                           required_only=args.required_only)
    print(f'{len(samples)} samples loaded', file=sys.stderr)

    probes = {}
    if args.structural_variants:
        print(f'loading structural variants from {args.structural_variants} max maf = {args.max_sv_maf}',
              file=sys.stderr)
        probes, structural_variants = simulation_noise.load_structural_variants(args.structural_variants,
                                                                                max_maf=args.max_sv_maf)
        print(f'{len(structural_variants)} svs loaded', file=sys.stderr)
        print(f'{len(probes)} probes after loading svs', file=sys.stderr)
        
    if args.confounding_snps:
        print(f'loading snps from {args.confounding_snps} max maf = {args.max_snp_maf}, max_distance = {args.max_distance}', file=sys.stderr)
        probes, confounding_snps = simulation_noise.load_confounding_snps(args.confounding_snps, 
                                                                          probes=probes, 
                                                                          max_distance=args.max_distance,
                                                                          max_maf=args.max_snp_maf)
        print(f'{len(confounding_snps)} snps loaded', file=sys.stderr)
        print(f'{len(probes)} probes after loading snps', file=sys.stderr)
        
    if args.probe_stats:
        probes = simulation_noise.load_probe_stats(args.probe_stats, probes=probes)
        print(f'{len(probes)} probes after loading stats', file=sys.stderr)

    out_file = sys.stdout
    if args.output:
        out_file = matrix_utils.open_output_file(args.output)

    include_meta = False
    header = matrix_utils.get_header_from_file(args.input,
                                               include_meta,
                                               required=required, 
                                               start_idx=args.probe_start_idx,
                                               required_only=args.required_only)
    
    
    print(header, file=out_file)
    
  
    counter = make_combinations(samples, 
                                args.choose, 
                                args.with_replacement,
                                out_file,
                                args.max_individuals,
                                args.tissue,
                                args.stage,
                                probes)
    out_file.close()
    print(f'{counter} simulated samples written', file=sys.stderr)
    print('completed', file=sys.stderr)
 

if __name__ == "__main__":
    main()    
