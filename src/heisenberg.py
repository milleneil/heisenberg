#! /usr/bin/env python3

import sys

import calculate_probe_stats as calc_probe
import simulate_samples as sim_sample
import simulate_methyl_mixture as sim_mix
import extract_sra_probe_vals as extract_probe
import extract_sra_sample_metadata as extract_meta
import invert_sra_table as invert_table
import subset_probes_samples as subset
import combine_sra_projects as combine
import print_cell as cell

def usage():
    print(
    '''
usage: heisenberg <module> <options>
    
To see help text for individual module, run:

heisenberg.py <module>    
    
  ---- simulation --------------------------------------------------------------
  simulate               simulate samples
  mix                    simulate cfDNA mixtures of two tissue types

  ---- preparation -------------------------------------------------------------
  stats                  read probe file and gather descriptive statistics 
  extract_sra_probe      extract SRA probe values from series_matrix.txt files
  invert                 turn 'long' file (probes as rows) into 'wide' (probes
                         as columns)

  ---- utility -----------------------------------------------------------------
  extract_sra_meta       extract SRA sample metadata from series_matrix.txt
  subset                 extract probe or sample subset from master file
  combine                safely combine multiple files into one
  cell                   print value of master file cell[x,y]                       
    '''
    )


modules = { 
    'stats' : calc_probe,
    'simulate' : sim_sample,
    'mix' : sim_mix,
    'extract_sra_probe' : extract_probe,
    'invert' : invert_table,
    'extract_sra_meta' : extract_meta,
    'subset' : subset,
    'combine' : combine,
    'cell' : cell }

if len(sys.argv) < 2 or sys.argv[1] not in modules:
    usage()
    sys.exit(1)

module = sys.argv.pop(1)
func = modules[module]
func.main()