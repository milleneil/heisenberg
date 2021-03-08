
import matrix_utils
import sys

# column index of first probe value if demographic info not in file
PROBE_START_IDX = 4


class MethylSample:
    def __init__(self, probe_label_idxs=None, start_idx=4, line=None, required=None, missing=-1.0, required_only=False):
        """ 
        initialize sample from line if provided - use probe_label_idxs map
        to define which probe labels are in which columns; start_idx identifies
        first probe column; required = set of probe labels to include no
        matter what, sub in missing val for invalid/missing probes
        """
        
        # record probe vals indexed by label
        self.probe_vals = {}
        
        if line == None:
            self.case = None
            self.sample = None
            self.biospecimen = None
            self.tissue = None
        else:
            # need probe label indexes set if initializing by line
            if probe_label_idxs == None:
                raise Exception("label indexes required when init by line")
            fields = line.rstrip().split('\t')
            self.case = fields[0]
            self.sample = fields[1]
            self.biospecimen = fields[2]
            self.tissue = fields[3]
            
            if start_idx > 4:
                self.gender = fields[4]
                self.age = float(fields[5])
                self.age_group = fields[6]
                self.stage = fields[7]
            
            for i in range(start_idx, len(fields)):
                # default missing vals to -1 to differentiate between true zero
                label = probe_label_idxs[i]

                # skip unneeded probes if limiting to required only
                if required != None and required_only and label not in required:
                    continue
                    
                try:
                    val = float(fields[i])
                    
                except ValueError:
                    val = missing
                self.probe_vals[label] = val
        
            # set placeholders for any required probes not in file
            if required != None:
                for r in required:
                    if r not in self.probe_vals:
                        self.probe_vals[r] = missing
          
    def __str__(self): 
        print_vals = [self.case, self.sample, self.biospecimen, self.tissue]
        
        # include metadata if set
        if hasattr(self, 'gender'):
            print_vals.extend([self.gender, 
                               str(self.age), 
                               get_age_group(self.age), 
                               self.stage])
        
        # write probes in alphabetical order
        for p in sorted(self.probe_vals.keys()):
            #print_vals.append(self.sample + "_" + p)
            print_vals.append('{0:0.7f}'.format(self.probe_vals[p]))
        return '\t'.join(print_vals)
        
def load_probe_label_indexes(file, start_idx=4):        
    """
    read header line of file to get dict of probe labels referencing column
    index for that probe
    """
    f = matrix_utils.open_file(file)
    
    fields = f.readline().rstrip().split('\t')
    probe_label_idxs = {}
        
    # record probe label indexes so we can refer to probes by name
    for i in range(start_idx, len(fields)):
        probe_label_idxs[i] = fields[i]
        
    f.close()
    
    return probe_label_idxs

def load_file_as_samples(file, required=None, start_idx=4, required_only=False):
    """ 
    load wide methyl file and return as collection of MethylSamples
    
    ensure that any probes specified in required list are included in
    sample with placeholder val if missing 
    """

    # save with sample identifier as key referencing sample obj
    samples = {}

    # get probe label indexes so we can refer to probes by name
    probe_label_idxs = load_probe_label_indexes(file, start_idx=start_idx)

    # open file and read past header line
    f = matrix_utils.open_file(file)
    f.readline()
    
    # convert all others into methyl sample
    counter = 0
    for line in f:
        if line.strip() == "":
            continue
        counter += 1
        if counter % 100 == 0:
            print(f'{counter} lines loaded', file=sys.stderr)
        sample = MethylSample(probe_label_idxs, start_idx=start_idx, line=line, required=required, required_only=required_only)
        samples[sample.sample] = sample
    f.close()    
    return samples       

def get_age_group(age):
    """function to encapsulate assignment of age to our pre-defined age groups"""
    group = None
    if age < 35:
        group = "0-34"
    elif age < 55:
        group = "35-54"
    elif age < 65:
        group = "55-64"
    elif age < 75:
        group = "65-74"
    else:
        group = "75+" 

    return group

def load_sra_demographics(file, samples):
    """load demographic data for SRA samples"""

    # read file and set age and gender as values in sample
    with open(file, 'r') as f:
        # 0 project
        # 1 sample
        # 2 age
        # 3 age_group - not used, calculate using current logic
        # 4 gender
        
        # keep track of samples found to make sure we got metadata for everyone
        found = set()
        
        for line in f:
            # read header to figure out if demographic info is included in
            # file or not
            fields = line.rstrip().split('\t')
            
            if line.startswith('project'):
                continue
            if fields[1] in samples:
                sample = samples[fields[1]]
                sample.age = float(fields[2])
                #age_group = fields[3] -- ignore, determine dynamically
                sample.gender = fields[4]
                found.add(fields[1])
                
        for sample_label, sample in samples.items():
            if sample_label not in found:
                print(f'metadata not found for sample: {sample_label}', file=sys.stderr)
                
def load_tcga_demographics(file):
    """
    load demographic data for tumor (TCGA-PAAD) samples - expected to be 
    clinical.project/clinical.tsv
    """
    
    # return as dict with sample id as key referencing tuple of age, group,
    # gender and tumor stage
    metadata = {}
    
    with open(file, 'r') as f:
        for line in f:
            # read past header
            if line.startswith('case'):
                continue

            fields = line.rstrip().split('\t')
            submitter_id = fields[1]
            gender = fields[3]
            tumor_stage = fields[11]
                        
            age_in_days = int(fields[12])
            age = age_in_days / 365
            age_group = get_age_group(age)
            
            metadata[submitter_id] = (age,age_group,gender,tumor_stage)
    
    return metadata

def set_tcga_demographics(samples,meta):
    """
    read demographic data from meta dict and set as member variables on
    tcga sample instance
    """
    for sample_label, sample in samples.items():
        age,_,gender,stage = meta[sample.case]
        sample.age = int(age)
        sample.gender = gender
        sample.stage = stage


def match_by_demographics(samples, age, gender):
    """get demographically matched samples for age group and gender"""
    matched = {}
    
    age_group = get_age_group(age)
    for sample_id, sample in samples.items(): 
        s_age_group = get_age_group(sample.age)   
        if s_age_group == age_group and gender == sample.gender:
            matched[sample_id] = sample
    return matched

def filter_by_demographics(samples, min_age):
    """
    filter sample set to remove samples with no demographic metadata and
    those younger than minimum age
    """
    filtered = {}
    for sample_id, sample in samples.items():
        if hasattr(sample, 'age') and hasattr(sample, 'gender') and sample.age >= min_age:
            filtered[sample_id] = sample
        else:
            print(f'filtering sample by metadata {sample_id}', file=sys.stderr)
            
    return filtered

def write_samples_to_file(samples, file_name, header):
    """
    write submitted samples to a tmp file - this
    is done to persist subsets of samples temporarily so that they can
    be processed in multiple threads without having to copy the entire dataset
    between processes (as python multiprocessing does)
    """
    f = matrix_utils.open_output_file(file_name)
    print(header, file=f)
    for s in samples:
        if not hasattr(s, 'stage'):
            s.stage = 'tmp_stage'
        print(s, file=f)
    f.close()
    
    
    