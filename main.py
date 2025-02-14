#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this calls for all the script that extract features from the sequences
and calls for machine learning models for the classification
"""
import yaml
from sys import argv
from os.path import exists
import subprocess as sp

from general_functions import dir_maker
import blast
import diamond_filterer
import orphan_selector
import diamond_feature_extractor
import seq_feature_extractor
import machine_learning

def config_parser(config_path):
    """This function gets the input for the pipeline from the config file

    Args:
        config_path (str): path to the "config.yaml" file

    Returns:
        input (dict): input files and directories
        parameters (dict): parameters for tools in the pipeline
        output_directories (dict): the names of the output directories
        negative_directories (dict): all the output directories for the negative set
    """
    with open(config_path) as ymlfile:
        cfg = yaml.full_load(ymlfile)
    
    input = cfg['input']
    parameters = cfg['parameters']
    output_directories = cfg['output_directories']
    return input, parameters, output_directories

def sim_seq_to_known_prot(known_prot, query, out_dir, n_threads, blast_type = 'blastp', eval = [0.001]):
    """This function maps the simulated sequences back to the known proteins using the "blast.py" script

    Args:
        known_prot (str): path to known protein fasta file
        query (str): path to fasta file containing the query sequences
        out_dir (str): the directory in which the files will be written
        n_threads (int): the number of threads that should be used
        blast_type (str, optional): the type of blast that should be used. Defaults to 'blastp'.
        eval (list, optional): the cut-off for the evalue. Defaults to [0.001].

    Returns:
        str: path to the output directory containing the blast results
    """
    #check if the directory with the highest evalue already exists in this case the 
    #assumption is made that the diamond has already been executed
    if exists(f'{out_dir}_eval{max(eval)}'):
        print(f'{out_dir}_eval{max(eval)} already exists, therefore no diamond run will be executed')
    else:   
        #make a directory to put the results into
        dir_maker(out_dir)
        #Use diamond to align the simulated sequences against the entire protein dataset
        out_file = f"{out_dir}/diamond.m12"
        #Run the blast using the blast script
        blast.main(known_prot, query, blast_type, out_file, n_threads, max(eval))
        #filter the table further into different evalue cut-offs
        diamond_filterer.main([out_file, eval])
    return out_dir

def main():
    inp, prmt, out_dir = config_parser(argv[1])
    #map the simulated proteins back to database of known proteins
    sim_seq_to_known_prot(inp['in_fasta'], inp['orphan_db'], out_dir['out_diamond'], inp['n_threads'], eval=prmt['evalues'])
    #select the orphans in a certain range of evalues
    diamond_table, diamond_seq = orphan_selector.main(out_dir['out_diamond'], inp['orphan_db'], [prmt['evalues'][0], prmt['evalues'][-1]], inp['n_threads'], out_dir['out_orphan'])
    #extract features from the diamond table and sequences for the machine learning model
    diamond_features = diamond_feature_extractor.main(inp['orphan_db'], diamond_table, diamond_seq, out_dir['out_machine_learning'])
    seq_features = seq_feature_extractor.main(out_dir['out_machine_learning'], diamond_features, inp['orphan_db'])
    classification = machine_learning.main(diamond_features, seq_features, inp['ml_model'])
    with open('test.txt', 'w') as f:
        f.write(classification)

if __name__ == '__main__':
    main()