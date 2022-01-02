#! /usr/bin/env python3

import numpy as np
import pandas as pd
import gzip
import logging


def load_labels(input_file, tissue_idx=3, numeric=False, as_is=False):
    """ 
    Read 'wide' input file and extract vector of tissue labels. 
    By default, assume tissue type is column idx 3.  Code with:
        0 = Normal
        1 = Tumor

    If as_is flag is True, return textual labels as they appear
    else resolve into two simple classes:
        Normal
        Tumor (where Tumor = everything not Normal) 
    """
    y = []
    f = open_file(input_file)
    first = True

    label_nums = get_labels_number_map()

    for line in f:
        if first:
            first = False
            continue
        fields = line.rstrip().split('\t')
        if numeric:
            if 'Normal' in fields[tissue_idx] or 'normal' in fields[tissue_idx]:
                y.append(label_nums['Normal'])
            else:
                y.append(label_nums['Tumor'])
        elif as_is:
            y.append(fields[tissue_idx])
        elif 'Normal' in fields[tissue_idx] or 'normal' in fields[tissue_idx]:
            y.append('Normal')
        else:
            y.append('Tumor')
    f.close()
    return np.array(y)


def num_label_for(label):
    if 'Normal' in label or 'normal' in label:
        return 1
    else:
        return 0


def get_number_labels_map():
    """
    Get dictionary that defines our number -> label encoding for easy lookup 
    """
    return {0: 'Tumor', 1: 'Normal'}


def get_labels_number_map():
    """
    Get dictionary that defines our label -> number encoding for easy lookup 
    """
    return {'Tumor': 0, 'Normal': 1}


def load_numeric_labels(input_file, tissue_idx=3):
    """
    Read 'wide' input file and extract vector of tissue type labels
    as numbers. Convenience method wrapping load_labels
    """
    return load_labels(input_file, tissue_idx, numeric=True)


def load_file_to_matrix(input_file, cols=None, probe_start=4):
    """ Load probe alpha values into matrix """
    # skip first three columns that have sample metadata and just keep
    # numeric vals
    f = open_file(input_file)
    # peek at the first line of file to get probe labels
    # from header and figure out how many columns there are
    # so we can skip itemizing them
    fields = f.readline().rstrip().split('\t')

    # if columns of interest not defined, read all columns beyond
    # sample/biospecimen/tissue info
    if cols is None:
        num_cols = len(fields)
        labels = fields[probe_start:]
        f.seek(0)
        cols = range(probe_start, num_cols)
    else:
        # extract out labels for probes we want
        labels = []
        for idx in cols:
            labels.append(fields[idx])
    f.close()
    matrix = np.genfromtxt(input_file,
                           delimiter='\t',
                           skip_header=1,
                           usecols=cols,
                           filling_values=0)
    return matrix, labels


def get_column_labels(input_file, cols=None, probe_start=3):
    # skip first three columns that have sample metadata and just keep
    # numeric vals
    labels = None
    f = open_file(input_file)
    # peek at the first line of file to get probe labels
    # from header and figure out how many columns there are
    # so we can skip itemizing them
    fields = f.readline().rstrip().split('\t')
    f.close()

    # if columns of interest not defined, read all columns beyond
    # sample/biospecimen/tissue info
    if cols == None:
        labels = fields[probe_start:]
    else:
        # extract out labels for probes we want
        labels = []
        for idx in cols:
            labels.append(fields[idx])

    return labels


def open_file(input_file):
    """ convenience method to open regular and gzipped files """
    if input_file.endswith('.gz'):
        f = gzip.open(input_file, 'rt')
    else:
        f = open(input_file, 'r')

    return f


def open_output_file(input_file):
    """ convenience method to open regular and gzipped files """
    if input_file.endswith('.gz'):
        f = gzip.open(input_file, 'wt')
    else:
        f = open(input_file, 'w')

    return f


def get_header_from_file(in_file, has_metadata, required=None, start_idx=4,
                         required_only=False):
    """
        get header from input file and print to output file -  
        insert metadata and required cols if provided - probe labels will
        be printed in sorted order regardless of input order
        """
    f = open_file(in_file)
    line = f.readline()
    f.close()
    columns = line.rstrip().split('\t')

    # make header - include initial cols up to first probe no matter what
    to_print = columns[:start_idx]
    if has_metadata:
        # add on metadata if provided
        to_print.extend(['gender', 'age', 'age_group', 'tumor_stage'])

    # add on probe columns - sort here to make sure we line up with
    # sorted vals in output samples
    if required != None and required_only:
        to_print.extend(sorted(required))
    else:
        probe_cols = set(columns[start_idx:])
        if required != None:
            probe_cols.update(required)
        to_print.extend(sorted(probe_cols))
    return '\t'.join(to_print)


def indexes_for_labels(input_file, labels):
    """ 
    read header of input file and find column indexes for submitted
    column labels.  Returned list will contain indexes sorted 
    numerically regardless of order they appear in input list 
    """
    f = open_file(input_file)
    fields = f.readline().rstrip().split('\t')
    f.close()

    col_idxs = []

    # keep track of found labels 
    found = set()
    for i in range(0, len(fields)):
        if fields[i] in labels:
            col_idxs.append(i)
            found.add(fields[i])

    # complain loudly if any labels weren't found in file
    for label in labels:
        if label not in found:
            raise Exception("label not in file: " + label)

    return sorted(col_idxs)


def load_probe_list(file):
    """ load list of probe names from file into list """
    f = open_file(file)
    probes = set()
    for line in f:
        fields = line.rstrip().split('\t')
        probes.add(fields[0])
    f.close()
    return probes

def load_as_dataframe(input_files, probes=None):
    """
    Load one or more matrix files into pandas dataframe, optionally restricting
    to list of probes
    """
    files = []

    # make sure to include sample and tissue type in addition to probe names
    if probes is not None:
        probes.update(['sample', 'tissue'])
    for f in input_files:
        logging.info(f'loading input file {f}...')

        df = pd.read_csv(f, sep='\t', usecols=probes)
        files.append(df)
        logging.info(f'loaded df {df.shape} from {f}...')

    if probes is not None:
        [probes.discard(f) for f in ['sample', 'tissue']]

        return pd.concat(files)