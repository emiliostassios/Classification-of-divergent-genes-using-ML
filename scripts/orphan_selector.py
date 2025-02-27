#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script selects all the proteins that did not get a diamond hit
with an evalue of the first evalue given in the range and did get a hit with the second
evalue given in the evalue range
"""

import os
from pathlib import Path
import re
import subprocess as sp
from general_functions import dir_maker
from Bio.SeqIO.FastaIO import FastaIterator, as_fasta
from concurrent.futures import ThreadPoolExecutor

def query_selector(diamond_table):
    """This function gets the query sequence which is the first column of the table

    Args:
        diamond_table (str): path to the diamond table

    Returns:
        list: query ids from the first column
    """
    query_ids = []
    with open(diamond_table) as f:
        for line in f:
            query_ids.append(line.split('\t')[0])
    unique_query_ids = list(set(query_ids))
    return unique_query_ids

def diamond_table_comparer(low_eval_table, high_eval_table):
    """This function compares two diamond tables and select the query ids
    that are present in the table with a high evalue cutoff and are not present
    in the table with a low evalue cutoff
    
    Args:
        low_eval_table (str): path to diamond table with low evalue cutoff
        high_eval_table (str): path to diamond table with high evalue cutoff

    Raises:
        ValueError: when the evalue of the high_table is lower than the evalue of the
        low_table

    Returns:
        dict: ids that are present in the high table and not in the low table
    """
    #get the evalues from the paths
    low_eval = float('.'.join(re.findall(r'\d+', re.split(r'/', low_eval_table)[0])))
    high_eval = float(re.findall(r'\d+', re.split(r'/', high_eval_table)[0])[0])
    #get query ids from both files
    if high_eval < low_eval:
        raise ValueError(f"The evalues are not in ascending order.\n \
           This raised the problem that the orphan sequences are selected\n \
           which don't get a hit with {low_eval_table} and do get hit with an\n \
           eval of {high_eval_table}")
    else:
        low_eval_ids = set(query_selector(low_eval_table))
        high_eval_ids = set(query_selector(high_eval_table))
        #get the ids that are present in high_eval_ids and not present in low_eval_ids
        filtered_ids = {id:None for id in high_eval_ids if id not in low_eval_ids}
        return filtered_ids

def seq_selector(in_seq, out_seq, ids):
    """This function gets the sequences corresponding to the ids given

    Args:
        in_seq (str): the input fasta file from which sequences should be selected
        out_seq (str): path to the output file
        ids (dict): the ids of the sequences that should be selected
    """
    if os.path.exists(out_seq):
        print(f'{out_seq} already exists, moving on...')
    else:
        with open(in_seq) as handle:
            for record in FastaIterator(handle):
                if record.id in ids:
                    with open(out_seq, 'a') as handle_out:
                        handle_out.write(as_fasta(record))
                    
        print(f'{out_seq} has been made')

def line_query_selector(diamond_table, query_ids, out_table):
    """This function selects all the lines from a diamond table that has the
    query ids

    Args:
        diamond_table (str): path to diamond table
        query_ids (list): query ids that should be selected from the diamond table
        out_table (str): path to the output table
    """
    if os.path.exists(out_table):
        print(f'{out_table} already exists, moving on...')
    else:
        #get all lines that have the query id
        lines = []
        with open(diamond_table) as f:
            for line in f:
                query_id = line.split('\t')[0]
                if query_id in query_ids:
                    lines.append(line)
        #write the lines to a file
        with open(out_table, 'w') as f:
            for line in lines:
                f.write(line)
        print(f'{out_table} has been made')

def table_filterer(args):
    """This function removes all the diamond hits from the fasta 
    and write it to another file

    Args:
        args (list): containg 5 elements
            element 1 (str): path to the diamond table with low eval cutoff
            element 2 (str): path to the diamond table with high eval cutoff
            element 3 (str): path to the output file for sequences
            element 4 (str): path to output file for diamond results
            element 5 (str): path to the input simulated sequences
    """
    sim_fa = args[0]
    low_table = args[1]
    high_table = args[2]
    out_seq = args[3]
    out_table = args[4]
    if os.path.exists(out_seq) and os.path.exists(out_table):
        print(f'"{out_seq}" and "{out_table} already exists, moving on...')
    else:
        #get ids that are present in the high_table and not in the low_table
        ids = diamond_table_comparer(low_table, high_table)
        print(ids)
        #write the sequences that are not in the ids to a seperate file
        seq_selector(sim_fa, out_seq, ids)
        #filter the ids from the diamond table
        line_query_selector(high_table, ids, out_table)

def main(diamond_dir, orphan_seq, eval_range, threads, out_dir):
    """This function selects the sequences that don't have a diamond hit 
    at eval1 and do have a hit at eval2

    Args:
        diamond_dir (str): path to directory containing the diamond results
        orphan_seq (str): path to file containing the orphan sequences
        eval_range (list): two evalues, sequences will be selected that do
        have a hit a the first value and don't have a hit at the second value
        threads (int): the number of threads you want to use
        out_dir (str): path to directory containing the sequences with no hit

    Returns:
        str: path to diamond and the corresponding sequences
    """
    args_table_filterer = []

    main_dir_eval1 = f'{diamond_dir}_eval{eval_range[0]}'
    main_dir_eval2 = f'{diamond_dir}_eval{eval_range[1]}'
    for file in os.listdir(f'{main_dir_eval1}'):
        path_file1 = f'{main_dir_eval1}/{file}'
        path_file2 = f'{main_dir_eval2}/{file}'
        #make directory to put the results into
        dir_maker(out_dir)
        out_file_seq = f'{out_dir}/orphan.fa'
        out_file_table = f'{out_dir}/diamond_orphan.m12'
        args_table_filterer.append([orphan_seq, path_file1, path_file2, out_file_seq, out_file_table])
    #Remove all the diamond hits from the fasta and write it to another file
    with ThreadPoolExecutor(threads) as thread:
            thread.map(table_filterer, args_table_filterer)
    return out_file_table, out_file_seq

if __name__ == '__main__':
    main()