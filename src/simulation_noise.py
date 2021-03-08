"""
Collection of utilities for adding noise/confounding simulated methyl values
"""

import matrix_utils
import random
import sys

class Probe:
    """
    class to represent methylation probe
    """
    
    def __init__(self, label, **kwargs):
        self.label = label
        self.chrom = kwargs.get('chrom')
        self.start = kwargs.get('start')
        self.stop = kwargs.get('stop')
        self.mean = kwargs.get('mean')
        self.stdev = kwargs.get('stdev')
        self.minimum = kwargs.get('minimum')
        self.maximum = kwargs.get('maximum')
        self.svs = []
        self.snps = []
        
    def __str__(self):
        return self.label + ' pos[' + str(self.chrom) + ':' + str(self.start) \
                    + '-' + str(self.stop) \
                    + '] stats[min=' + str(self.minimum) + ',max=' + str(self.maximum) \
                    + ',mean=' + str(self.mean) + ',stdev=' + str(self.stdev)


    @property
    def start(self):
        return self.__start
    
    @start.setter
    def start(self, start):
        self.__start = safe_int(start)
    
    @property                           
    def stop(self):       
        return self.__stop
     
    @stop.setter
    def stop(self, stop):
        self.__stop = safe_int(stop)

    @property                           
    def minimum(self):       
        return self.__minimum
     
    @minimum.setter
    def minimum(self, minimum):
        self.__minimum = safe_float(minimum)

    @property                           
    def maximum(self):       
        return self.__maximum
     
    @maximum.setter
    def maximum(self, maximum):
        self.__maximum = safe_float(maximum)

    @property                           
    def mean(self):
        return self.__mean
     
    @mean.setter
    def mean(self, mean):
        self.__mean = safe_float(mean)

    @property                           
    def stdev(self):
        return self.__stdev
     
    @stdev.setter
    def stdev(self, stdev):
        self.__stdev = safe_float(stdev)


class ProbeVariant:
    """
    class to represent a variant that confounds methylation probe
    """
    
    def __init__(self, label, **kwargs):
        self.label = label
        self.chrom = kwargs.get('chrom')
        self.start = safe_int(kwargs.get('start'))
        self.stop = safe_int(kwargs.get('stop'))
        self.type = kwargs.get('type')
        self.freq = safe_float(kwargs.get('freq'))
        self.hom_freq = safe_float(kwargs.get('hom_freq'))
        self.het_freq = safe_float(kwargs.get('het_freq'))
        self.probes = []
        
    def __str__(self):
        return self.label + ' pos[' + str(self.chrom) + ':' + str(self.start) \
                    + '-' + str(self.stop) + ', type=' + type

def safe_int(val):
    if val == None:
        return val
    else:
        return int(val)

def safe_float(val):
    if val == None:
        return val
    else:
        return float(val)


def load_probe_stats(file, probes = None):
    """ 
    read file of descriptive statistics for probes and return 
    dict of probe name referencing probe obj
    """
    
    # if probes hasn't already been initialized, do it now
    if probes == None:
        probes = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('probe'):
                continue
            label,f_min,f_max,f_mean,f_stdev = line.rstrip().split('\t')

            if label not in probes:
                probes[label] = Probe(label)
            probes[label].mean = f_mean
            probes[label].minimum = f_min
            probes[label].maximum = f_max
            probes[label].stdev = f_stdev
            
    return probes



def load_structural_variants(file, probes=None, max_maf=.3):
    """
    read file of structural variants and affected probes and return
    as dict of probe variant instances with label reference instance
    """
    svs = {}
    if probes == None:
        probes = {}
    
    f = matrix_utils.open_file(file)
    for line in f:
        # 0 probe
        # 1 chr
        # 2 probe start
        # 3 probe stop
        # 4 sv name
        # 5 sv type
        # 6 sv start
        # 7 sv stop
        # 8 sv allele frequency
        fields = line.rstrip().split('\t')
        probe_label = fields[0]
        label = fields[4]
        if label not in svs:
            svs[label] = ProbeVariant(label, 
                                      chrom = fields[1], 
                                      start = fields[6], 
                                      stop = fields[7],
                                      type = fields[5],
                                      freq = fields[8],
                                      hom_freq = fields[9],
                                      het_freq = fields[10])
        if probe_label not in probes:
            probes[probe_label] = Probe(probe_label,
                                        chrom = fields[1],
                                        start = fields[2],
                                        stop = fields[3])
        
        # exclude common variants
        if svs[label].freq <= max_maf:
            
            # set cross reference so each probe/variant pair knows about each other
            probes[probe_label].svs.append(svs[label])
            svs[label].probes.append(probe_label)
            
        
    f.close()
    
    return probes, svs

def load_confounding_snps(file, probes=None, max_distance=2, max_maf=.2):
    """
    read file of snps that may confound methyl probe values and return
    as probe id referencing snps
    """
    snps = {}
    if probes == None:
        probes = {}
    
    f = matrix_utils.open_file(file)
    for line in f:
        # 0 target id (probe)
        # 1 snp (; delimited)
        # 2 distance (from probe, ; delimited)
        # 3 minor allele frequency (; delimited)
        if line.startswith('TargetID'):
            continue
        probe_label, rsIds, distances, mafs = line.rstrip().split('\t')
        rsIds = rsIds.split(';')
        mafs = mafs.split(';')
        distances = distances.split(';')
        if probe_label not in probes:
            probes[probe_label] = Probe(probe_label)
        for i in range(0, len(rsIds)):
            # only keep snps within requested distance and that are rarer
            # than cutoff
            if int(distances[i]) <= max_distance and float(mafs[i]) <= max_maf:
                rsId = rsIds[i]
                if rsId not in snps:
                    snps[rsId] = ProbeVariant(rsId, freq = mafs[i])
                snps[rsId].probes.append(probe_label)
                probes[probe_label].snps.append(snps[rsId])
        
    f.close()
    
    return probes, snps


def get_probe_noise(probe):
    """ 
    get random amount of noise to apply to simulated values for submitted
    probe by selecting a random value with distribution defined by stdev
    for probe
    """
    noise = 0
    if probe.stdev != None:
        noise = random.uniform(-probe.stdev, probe.stdev) 
 
    return noise

def get_random_variants(variants, by_zygosity=False, homozygous=False):
    """ 
    identify random list of structural variants based on their MAFs by using 
    MAF as a probability that an individual has the variant. By default, use
    variant MAF - if requested, specifically get random hom vs. het variant
    based on zygosity specific allele frequency 
    
    This method obviously ignores any kind of linkage and assumes each variant 
    appears independently
    """
    
    chosen = set()
    
    for variant in variants:
        # randomly choose a number between 0 and 1 - if number is less than
        # variant allele frequency, keep it
        choice = random.random()
        if by_zygosity:
            compare_val = variant.hom_freq if homozygous else variant.het_freq
        else:
            compare_val = variant.freq
            
        if choice <= compare_val:
            chosen.add(variant)
            
    return chosen

def adjust_by_structural_variant(val, probe):
    """
    randomly test for sv that overlaps probe -- adjust value accordingly. For
    now, this is limited to dropping probe val to 0 if a hom. deletion is
    found
    """
    adjusted = val
    if len(probe.svs) > 0 :
        for sv in get_random_variants(probe.svs, by_zygosity=True, homozygous=True):
            if sv.type == 'DEL':
                print(f'zero out from del: {probe.label}, {sv.label}', file=sys.stderr)
                # set probe val to zero if we have a homozgyous deletion
                adjusted = 0
                break
                
            # how else to adjust vals? since beta val is an avg, we shouldn't
            # just double or halve value based on dup or het events (??
    return adjusted

def adjust_by_confounding_snp(val, probe):
    """
    randomly test for confounding snp, if found set probe val to 0
    """
    adjusted = val
    if len(probe.snps) > 0:
        snps = get_random_variants(probe.snps)
        if len(snps) > 0:
            snp = snps.pop()
            print(f'zero out from confounding snp: {snp.label} - {snp.freq}, probe: {probe.label}', file=sys.stderr)
            # if we have a snp, assume it confounds probe and set to zero
            adjusted = 0

    return adjusted

                
def random_adjust(val, probe):
    """ 
    perform all adjustments on probe value at once including stdev noise
    and simulation of sv and snv variants
    """
    
    # adjust by random noise (make zero min. possible)
    val = max(val + get_probe_noise(probe), 0)
    
    # drop probe val entirely if we have random structural variant or
    # confounding snp
    if val > 0:
        val = adjust_by_structural_variant(val, probe)
    
        # if we still have signal, see if rand snp would cause dropout
        if val > 0:
            val = adjust_by_confounding_snp(val, probe)
        
    return val

