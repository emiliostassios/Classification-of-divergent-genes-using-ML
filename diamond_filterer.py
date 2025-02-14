#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script is used to filter the output of diamond based on the e value
"""
from sys import argv
import re
import os
from general_functions import dir_maker
import subprocess as sp

def table_reader(in_table, cutoff):
    """This function reads in the table and removes lines that do not meet the evalue
    cutoff

    Args:
        in_table (str): path to the diamond table
        cutoff (int or float): the evalue cut-off

    Returns:
        list: all the lines that met the cutoff in a string format
    """
    fl_cutoff = float(cutoff)
    #Load in table and write the lines that meet the cut-off the the "lines" string
    lines = ''
    with open(in_table) as f:
        for line in f:
            eval = float(line.split(sep = '\t')[10])
            if eval < fl_cutoff:
                lines += line
    return lines

def main(args):
    """This function runs the main functions in the script

    Args:
        args (list): should contain 2 elements:
            element 1 (str): path to diamond output table (it does not matter which one)
            element 2 (int or float): represents the cut-off for the evalue

    Returns:
        str: the path to the output file. This is the same as the input path but '_eval..'
        (.. is the evalue cutoff) is added to the main directory in the path
    """
    diamond_table = args[0]
    eval_cutoffs = args[1]
    for cutoff in eval_cutoffs:
        #The highest cutoff is already used for the diamond
        if cutoff == max(eval_cutoffs):
            pass
        else:
            #make the directories to put the results into
            main_dir = diamond_table.split('/')[0]
            sub_dirs = re.sub(r'^.*?/', '/', diamond_table) 
            out_diamond = f'{main_dir}_eval{cutoff}{sub_dirs}'
            dirs = '/'.join(out_diamond.split('/')[:-1])
            #if the output file already exists move on to the next file
            if os.path.exists(out_diamond):
                print(f'The file "{out_diamond}" already exists, moving on...')
            else:
                out_lines = table_reader(diamond_table, cutoff)
                #Write the lines that met the cut-off to a seperate file
                dir_maker(dirs)
                
                with open(out_diamond, 'w') as f:
                    f.write(out_lines)
                    print(f'The file "{out_diamond}" was made, moving on...')
    #rename the original diamond directory which has the highest evalue
    main_dir = diamond_table.split('/')[0]
    sp.run(f'cp -r {main_dir} {main_dir}_eval{max(eval_cutoffs)}', shell=True)
    return out_diamond

if __name__ == '__main__':
    main()