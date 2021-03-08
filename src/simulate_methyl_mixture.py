#! /usr/bin/env python3

"""
Make simulated cell-free methylation values by mixing blood and tumor methylation
vals at requested concentration

Basic formula:

Bn = normal methylation val
Bt = tumor methylation val
Fn = fraction of normal (i.e. percent of total that should be normal)
Ft = fraction of tumor

Ba = Adjusted methyl value = (Fn x Bn) + (Ft x Bt)

with metadata file for both SRA and TCGA
    - divide samples by age group and gender
    - if requested, make normal set by matching


"""


import sys
import argparse
import matrix_utils
import methyl_sample_utils as m_utils
import random
import simulation_noise

def parse_args():
    parser = argparse.ArgumentParser(description='create simulated cell-free ' + 
                                     'DNA methylation probe values',
                                     prog='heisenberg mix')
    parser.add_argument('-n', '--normal', type=str, 
                        help='normal methylation values', required=True)

    parser.add_argument('-t', '--tumor', type=str, 
                        help='tumor methylation values', required=True)
    parser.add_argument('-a', '--all_by_all', action="store_true",
                        help='do all by all simulation')
 
    parser.add_argument('-o', '--output', type=str, 
                        help='output file to write [default=STDOUT]')
    
    parser.add_argument('-f', '--tumor_fraction', type=float, required=True, 
                        help='fraction of tumor signal in output (total=1)')

    parser.add_argument('-m', '--normal_metadata', type=str, 
                        help='metadata file with age and gender for normals')
    
    parser.add_argument('-u', '--tumor_metadata', type=str, 
                        help='metadata file with age and gender for tumors')
    
    parser.add_argument('-g', '--demographic', action="store_true", 
                        help='match samples by gender and age group')
    
    parser.add_argument('-e', '--min_age', type=int, 
                        help='exclude samples younger than this [default=35]',
                        default=35)

    parser.add_argument('-x', '--normal_probe_start_idx', type=int, default=4,
                        help='column index of first methyl probe in normal input [default=4]')
    
    parser.add_argument('-y', '--tumor_probe_start_idx', type=int, default=4,
                        help='column index of first methyl probe in tumor input [default=4]')
    
    parser.add_argument('-v', '--structural_variants', type=str,
                        help='structural variant probe overlap file')

    parser.add_argument('-c', '--confounding_snps', type=str,
                        help='list of SNPs that may confound methyl values')
    
    parser.add_argument('-s', '--probe_stats', type=str,
                        help='descriptive stats per probe for adding random noise')
    
    parser.add_argument('-p', '--probes',type=str,
                        help='only write values for these probes (skip all others)')

    
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


def adjust_vals(normal, tumor, normal_fraction, tumor_fraction, probes=None):
    adjusted = {}
    
    for probe_label in normal.keys():
        adjusted_normal = normal[probe_label] * normal_fraction
        adjusted_tumor = tumor[probe_label] * tumor_fraction
        combined_val = (adjusted_normal + adjusted_tumor)
        
        # simulate random noise if probes provided
        # randomly simulate sv
        # randomly simulate confounding snp
        # simulate random noise based off of stdev
        # adjust aggregate probe val accordingly
        if probes != None and probe_label in probes:
            probe = probes[probe_label]
            
            # negative values signal missing probe data in one or more
            # samples used to create simulated sample -- skip adjustment
            # in this case and print missing val (-1) for this probe
            if combined_val >= 0:
                combined_val = simulation_noise.random_adjust(combined_val, probe)
                adjusted[probe_label] = combined_val
            else:
                #FIXME: parameterize missing val
                adjusted[probe_label] = -1
        else:
            adjusted[probe_label] = combined_val
        
        
    return adjusted


def combine(normal, tumor, normal_fraction, tumor_fraction, probes=None):
    combined = m_utils.MethylSample()        
    
    normal_fraction_str = '{0:0.3f}'.format(normal_fraction)
    tumor_fraction_str = '{0:0.3f}'.format(tumor_fraction)

    combined.case = 'ct-sim'
    combined.sample = normal.case + "_" + tumor.case
    combined.biospecimen = normal_fraction_str + ":" + tumor_fraction_str
    combined.tissue = 'tumor'
    print(f'creating simulated sample {combined.sample}', file=sys.stderr)
    combined.probe_vals = adjust_vals(normal.probe_vals, tumor.probe_vals,normal_fraction, tumor_fraction,probes=probes)
    
    if hasattr(tumor, 'gender'):
        combined.gender = tumor.gender
        combined.age = tumor.age
        combined.stage = tumor.stage

    return combined

def main():
    args = parse_args()
    
    out_file = sys.stdout
    if (args.output):
        out_file = matrix_utils.open_output_file(args.output)
        
    if args.all_by_all:
        print('simulating all x all', file=sys.stderr)
    
    probes = {}
    if args.structural_variants:
        print(f'loading structural variants from {args.structural_variants} max maf = {args.max_sv_maf}', 
              file=sys.stderr)

        probes, structural_variants = simulation_noise.load_structural_variants(args.structural_variants,
                                                                                max_maf=args.max_sv_maf)
        print(f'{len(structural_variants)} svs loaded', file=sys.stderr)
        print(f'{len(probes)} probes after loading svs', file=sys.stderr)
        
    if args.confounding_snps:
        print(f'loading snps from {args.confounding_snps} max maf = {args.max_snp_maf}, max_distance = {args.max_distance}', 
              file=sys.stderr)

        probes, confounding_snps = simulation_noise.load_confounding_snps(args.confounding_snps, 
                                                                          probes=probes,
                                                                          max_distance=args.max_distance,
                                                                          max_maf=args.max_snp_maf)
        print(f'{len(confounding_snps)} snps loaded', file=sys.stderr)
        print(f'{len(probes)} probes after loading snps', file=sys.stderr)
        
    if args.probe_stats:
        probes = simulation_noise.load_probe_stats(args.probe_stats, probes=probes)
        print(f'{len(probes)} probes after loading stats', file=sys.stderr)
    
    probe_subset = None
    if args.probes:
        print(f'restricting output to probe subset: {args.probes}', file=sys.stderr)
        probe_subset = matrix_utils.load_probe_list(args.probes)
    
    print(f"reading normal from {args.normal}", file=sys.stderr)
    all_normal_samples = m_utils.load_file_as_samples(args.normal, 
                                                      start_idx=args.normal_probe_start_idx,
                                                     required=probe_subset, 
                                                     required_only=(probe_subset != None))    
            
    tumor_metadata = None       
    if args.demographic:
        print(f"matching by demographics", file=sys.stderr)
        if not args.tumor_metadata:
            raise Exception('must provide metadata for demographic matching')
        
        # if normal metadata not supplied as a file, it must be included
        # in input file itself -- if not we should see errors downstream
        if args.normal_metadata:
            print(f'loading normal metadata from {args.normal_metadata}', file=sys.stderr)
            m_utils.load_sra_demographics(args.normal_metadata, all_normal_samples)
            all_normal_samples = m_utils.filter_by_demographics(all_normal_samples, args.min_age)

        for sample_id, sample in all_normal_samples.items():
            if not hasattr(sample, 'age'):
                raise Exception("no metadata for normal sample: " + sample_id)
        print(f'loading tumor metadata from {args.tumor_metadata}', file=sys.stderr)
        tumor_metadata = m_utils.load_tcga_demographics(args.tumor_metadata)

    
    print(f"{len(all_normal_samples)} normal lines read", file=sys.stderr)

    normal_fraction = 1 - args.tumor_fraction

    
    # get probe label indexes to enable referring to them by name
    probe_label_idxs = m_utils.load_probe_label_indexes(args.tumor, start_idx=args.tumor_probe_start_idx)
    
    header = matrix_utils.get_header_from_file(args.tumor, 
                                               (tumor_metadata != None),
                                               start_idx=args.tumor_probe_start_idx,
                                               required=probe_subset, 
                                               required_only=(probe_subset != None))
    print(header, file=out_file)

    
    # input files are rows = samples, cols = probe vals
    # open file and read past header
    f = matrix_utils.open_file(args.tumor)
    counter = 0
    for line in f:
        if line.startswith('case'):
            continue
        elif 'Metastatic' in line:
            continue
        counter += 1
        if counter % 10 == 0:
            print(f"{counter} tumor lines read", file=sys.stderr)
        
        tumor = m_utils.MethylSample(line=line, 
                                     probe_label_idxs=probe_label_idxs, 
                                     start_idx=args.tumor_probe_start_idx,
                                     required=probe_subset, 
                                     required_only=(probe_subset != None))
        normal_samples = all_normal_samples
        
        if tumor_metadata:
            
            age,_,gender,stage = tumor_metadata[tumor.case]
            tumor.age = int(age)
            tumor.gender = gender
            tumor.stage = stage
            
            normal_samples = m_utils.match_by_demographics(all_normal_samples,
                                                           tumor.age,
                                                           tumor.gender)
            num_samples = len(normal_samples)
            if num_samples == 0:
                print(f'No normal matches found for tumor sample: {tumor.age}|{tumor.gender}',
                       file=sys.stderr)
                continue
        
        # if desired, choose a single random sample to mix with 
        if args.all_by_all:
            for normal_label, normal in normal_samples.items():
                combined = combine(normal, tumor, normal_fraction, args.tumor_fraction, probes=probes)
                print(combined, file=out_file)
        else:
            normal_label = random.choice(list(normal_samples.keys()))
            normal = normal_samples[normal_label]
            
            combined = combine(normal, tumor, normal_fraction, args.tumor_fraction, probes=probes)
            print(combined, file=out_file)
     
    out_file.close()
    
    
    

if __name__ == "__main__":
    main()    
