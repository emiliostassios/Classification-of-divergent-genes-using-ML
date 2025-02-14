#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script extract features from raw sequences
"""
from general_functions import dir_maker
import subprocess as sp
from os.path import exists
from Bio.SeqIO.FastaIO import FastaIterator, as_fasta
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def diamond_parser(file):
    """This function parses the diamond feature tables

    Args:
        file (str): paths to diamond table

    Returns:
        list: with query ids (first column)
    """
    query_ids = []
    with open(file) as f:
        next(f)
        for line in f:
            query_ids.append(line.split('\t')[0])
    return query_ids

def fasta_filterer(fasta, ids, prefix = None):
    """This function selects sequences that match the id list

    Args:
        fasta (str): path to fasta file
        ids (list): the record ids that should be selected
        prefix(str): The prefix that could be appended to the record id.
        Defaults to None

    Returns:
        list: all the selected records
    """
    #get all the records that match the ids
    records = []
    with open(fasta) as handle:
        if prefix == None:
            for record in FastaIterator(handle):
                if record.id in ids:
                    records.append(record)
        else:
            for record in FastaIterator(handle):
                if record.id in ids:
                    record.id = f'{prefix}{record.id}'
                    records.append(record)
    return records

def ids_to_fasta(args):
    """This function selects all the sequences from the file_ids and
    writes it to a file

    Args:
        list: containing 5 elements:
            element 1 (list): ids that should also be the id of the fasta file
            element 2 (str): path to the fasta file
            element 3 (str): path to the output file
    """
    ids = args[0]
    in_fasta = args[1]
    out = args[2]
    #get all the output paths
    if exists(out) == False:
        #filter fasta files
        all_records = fasta_filterer(in_fasta, ids)
        print(f'records gathered for {out}')
        #make the directory to put the results into
        dir = '/'.join(out.split('/')[:-1])
        dir_maker(dir)
        #write the records to a file
        with open(out, 'w') as handle:
            for record in all_records:
                fa = as_fasta(record)
                handle.write(fa)

def UniRep(args):
    """This function runs the unirep_feature_extractor.py script.
    This script extracts features from raw sequences

    Args:
        list: containing 3 elements:
            element 1 (str): path to a fasta file from which features should be extracted
            element 2 (str/int): the batch size that should be used for unirep
            element 3 (str): path to where the features should be written
    """
    fa_file = args[0]
    batch_size = args[1]
    out_feature = args[2]
    if exists(out_feature):
        print(f'"{out_feature}" already exists, moving on...')
    else:
        #Note that the unirep conda env should be used here
        cmd = f'nice -16 conda run -n unirep python ../unirep_feature_extractor.py {fa_file} {batch_size} {out_feature} {False}'
        print(cmd)
        sp.run(cmd, shell = True)
        print(f'feature are extracted from {fa_file} \nand can be found at {out_feature}. Moving on...')

def main(ml_dir, diamond_table, sequences):
    """This function runs the functions required to get the desired sequences and 
    extract features from these sequences

    Args:
        ml_dir (str): path to the machine learning directory
        diamond_table (str): path to diamond table containing the diamond features
        sequences (str): path to the file with the sequences that are also in the diamond table
    
    Returns:
        str: path to the directory containing the feature tables
    """
    batch_size_unirep = 15
    #specify output
    tmp_dir = 'temp_out'
    out_fa = f'{tmp_dir}/pos_out.fa'
    out_features = f'{ml_dir}/seq_features'
    #get sequences corresponding to the Diamond feature ids
    query_ids = diamond_parser(diamond_table)
    ids_to_fasta([query_ids, sequences, out_fa])
    #get features from the sequences
    args_unirep = []
    dir_maker(out_features)
    out_file = f'{out_features}/seq_features.tsv'
    args_unirep.append([sequences, batch_size_unirep, out_file])
    #use Unirep to extract features from the sequences
    print("WARNING: Unirep is used in this script. Download Unirep from https://github.com/churchlab/UniRep.git if you haven't already.")
    print('Delete line 10 and replace it with "from UniRep.data_utils import aa_seq_to_int, int_to_aa, bucketbatchpad" in the "unirep.py" script.')
    #note that unirep uses many threads, which is why I would advice to use a maximum of 4 processes here.
    UniRep(args_unirep[0])
    return out_file

if __name__ == '__main__':
    main()