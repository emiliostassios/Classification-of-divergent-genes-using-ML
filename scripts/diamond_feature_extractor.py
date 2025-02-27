#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script extract features from diamond tables
"""

from orphan_selector import query_selector
from Bio.SeqIO.FastaIO import FastaIterator
import subprocess as sp
from os.path import exists
import portion as P
import numpy as np
from general_functions import dir_maker
from concurrent.futures import ThreadPoolExecutor

def fasta_parser(fasta):
    """This function parses a fasta file

    Args:
        fasta (str): path to the fasta file containing the query sequences

    Returns:
        dict: key: record id; value: record sequence
    """
    #get the query sequence
    records = {}
    with open(fasta) as handle:
        for record in FastaIterator(handle):
            records[record.id] = record.seq
    return records

def value_inside_interval(interval):
    """This function counts the number that is inside the interval

    Args:
        interval (list): a list of lists containing the intervals

    Returns:
        int: the counted number inside the interval
    """
    #make empty interval
    I = P.empty()
    #append the list of lists into the interval
    for x,y in interval:
        I = I | P.closed(x,y) 
    #count the value inside the interval
    n = 0
    for interval in I:
        n += (interval.upper - interval.lower)
    return n

def feature_avg_calculator(table, query_ids, records):
    """This function calculates the number of alignments and the 
    average percent identity, evalue, bit score and coverage for a given query.

    Args:
        table (str): path to a diamond table
        query_ids (list): query ids of interest
        records (dict): key: record id (that should correspond to a query id)
        value: the corresponding sequence

    Returns:
        dict: key: query id and value: list containing the alignment count, avg identity, 
        avg evalue, avg bit score, avg coverage, avg number of matches, highest bit score,
        the longest alignment length and the class of the query.
    """
    features = dict()
    for query in query_ids:
        alignment_count = 0
        percent_identity = 0
        eval = 0 
        bit_score = 0
        n_matches = 0
        evalues = []
        top_identity = []
        bit_score_list = []
        n_matches_list = []
        query_coverage = []
        with open(table) as f:
            for line in f:
                if line.startswith(query):
                    #sum all the values for each line that starts with the query
                    alignment_count += 1               
                    percent_identity += float(line.split('\t')[2])
                    eval += float(line.split('\t')[10])
                    bit_score += float(line.split('\t')[11])
                    n_matches += float(line.split('\t')[3])
                    evalues.append(float(line.split('\t')[10]))
                    top_identity.append(float(line.split('\t')[2]))
                    bit_score_list.append(float(line.split('\t')[11]))
                    n_matches_list.append(float(line.split('\t')[3]))
                    query_coverage.append([float(line.split('\t')[6]), float(line.split('\t')[7])])
            if alignment_count != 0:
                #Calculate averages
                avg_identity = percent_identity / alignment_count
                avg_eval = eval / alignment_count
                avg_bit_score = (bit_score / alignment_count) / len(records[query])
                #calculate average coverage
                avg_n_matches = n_matches / alignment_count
                avg_coverage = avg_n_matches / len(records[query])
                #calculate min and top scores
                min_eval = min(evalues)
                top_alignment = max(top_identity)
                top_bit_score = max(bit_score_list) / len(records[query])
                top_n_matches = max(n_matches_list) / len(records[query])
                #calculate query coverage
                counted_query_coverage = value_inside_interval(query_coverage) / len(records[query])
                #append all features to the features dictionary
                features[query] = [
                    alignment_count, avg_identity, avg_eval, avg_bit_score, avg_coverage, min_eval, top_alignment,
                    avg_n_matches, top_bit_score, top_n_matches, counted_query_coverage
                ]
    return features

def no_hit_adder(sim_fa, feature_dict):
    """This function appends keys to a dict that are not present.
    The value will be a list of NaN values of the length of the first
    value in the feature_dict

    Args:
        sim_fa (str): path to fasta file
        feature_dict (dict): a dict to which value from the fasta file should be added

    Returns:
        dict: dictionary with the appended values
    """
    #make row
    for key,value in feature_dict.items():
        value_len = len(feature_dict[key]) - 1
        break
    row = [np.nan for i in range(value_len)]
    row.insert(0, 0)
    #get all the records from the simulated sequences
    records = fasta_parser(sim_fa)
    for id, seq in records.items():
        feature_dict.setdefault(id, row)
    return feature_dict         

def tsv_format_maker(header, features):
    """this function converts the features to tsv format

    Args:
        header (list): the header of the table
        features (dict): key: first column; value: other columns

    Returns:
        str: tsv formatted string
    """
    #convert the dictionaries to a tsv formatted list10
    values = ['{}\t{}\n'.format(key, "\t".join(map(str, value))) for key, value in features.items()]
    #merge everything into one big tsv formatted list
    tsv_format = header
    tsv_format.extend(values)
    return tsv_format

def diamond_feature_extractor(args):
    """This function extracts features from a diamond table for each
    given query id. For each query it calculates the number of alignments, 
    the average identity, average evalue, average bit score and the query coverage

    Args:
        args (list): containg 5 elements
            element 1 (str): path to the diamond table
            element 2 (str): path to the fasta file containing the query sequences
            element 3 (str): path to simulated fasta file
            element 4 (str): path to the output file
    """
    table = args[0]
    fasta = args[1]
    sim_fa = args[2]
    out = args[3]
    if exists(out):
        print(f'The file "{out}" already exists, moving on...')
    else:
        #Get all unique query ids
        query_ids = query_selector(table)
        #get sequences
        sequences = fasta_parser(fasta)
        #get features
        features_temp = feature_avg_calculator(table, query_ids, sequences)
        #append sequences that did not get a hit
        features = no_hit_adder(sim_fa, features_temp)
        #convert features to tsv format
        header1 = 'query\talignment_count\tavg_identity\tavg_eval\tavg_bit_score\tavg_coverage\tmin_eval\thighest_pident\t'
        header2 = 'avg_alignment_length\thighest_bit_score\thighest_alignment_length\tquery_coverage\n'
        header = [header1 + header2]
        tsv_format = tsv_format_maker(header, features)
        #make dir to put the table into
        dir = '/'.join(out.split('/')[:-1])
        dir_maker(dir)
        #write the list to a file
        with open(out, 'a') as f:
            for line in tsv_format:
                f.write(line)
        print(f'The file "{out}" has been made, moving on...')

def main(orphan_seq, diamond_table, diamond_seq, out_dir):
    """This function extracts features from a diamond table.
    Note that the diamond table should have the default table format:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    Args:
        orphan_seq (str): path to file containing the proteins sequences
        diamond_table (str): path to the file containing the diamond table
        diamond_seq (str): path to the file containing the sequences corresponding to the diamond table
        out_dir (str): path to where the output directory should be written
    Returns:
        str: path to the output directory of the features tables
    """
    out_features = f'{out_dir}/diamond_features'

    #get input for the "diamond_feature_extractor" function
    out = f'{out_features}/diamond_features.tsv'
    #add the paths to the list
    feature_extraction_args = [diamond_table, diamond_seq, orphan_seq, out]
    #extract features for each query id
    diamond_feature_extractor(feature_extraction_args)
    #Get the table with the least amount of rows and select this number of rows for all tables
    return out

if __name__ == '__main__':
    main()