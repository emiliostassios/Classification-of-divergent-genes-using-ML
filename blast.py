#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script runs diamond to align sequences agains a database
"""

from sys import argv
import subprocess as sp
import os
from pathlib import Path

def diamond_makedb(db):
    """This function makes a diamond db

    Args:
        db (str): path to fasta file from which a diamond database should be created

    Returns:
        str: path to the diamond db
    """

    db_name = 'diamond_' + Path(db).stem
    diamond_db = db_name + '.dmnd'
    if os.path.exists(diamond_db):
        print('Diamond database already exists, moving on...')
        return diamond_db
    elif db.endswith('.gz'):
        command = f'gunzip -c {db} | diamond makedb --in - -d {db_name}'
        sp.run(command, shell = True)
        return diamond_db
    else:
        command = f'diamond makedb --in {db} -d {db_name}'
        sp.run(command, shell = True)
        return diamond_db

def diamond_runner(db, query, blast_type, out, eval, threads = 20):
    """This function runs diamond

    Args:
        db (str): path to a diamond database
        query (str): path to the query for run
        blast_type (str): the blast type that should be run. The options are blastp or blastx
        out (str): path to the output of the run
        eval (int): cut-off for the e-value
        threads (int, optional): the number of threads that should used to run diamond. Defaults to 20.

    Returns:
        str: the path to where the output of the alignment is written
    """

    if os.path.exists(out):
        print("output file already exists, moving on...")
        return out
    command = f'nice -5 diamond {blast_type} -d {db} -q {query} --ultra-sensitive -e {eval} -o {out} -p {threads}'
    sp.run(command, shell = True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    print(f'Diamond run finished, file can be found at "{out}"')
    return out

def main(db, query, blast_type, out_file, n_threads, eval):
    """This function runs the main functions of this script

    Args:
        db (str): path to the database
        query (str): path to the query
        blast_type (str): the type of blast that needs to be run
        out_file (str): the path to the output for diamond
        n_threads (int): the number of threads that should be used
        eval (int or float): cut-off for the e-value.

    Returns:
        str: the path to the output for diamond
    """

    #Diamond only supports blastp and blastx
    if blast_type in ['blastp', 'blastx']:
        #make diamond database
        diamond_db = diamond_makedb(db)
        #run diamond
        alignment = diamond_runner(diamond_db, query, blast_type, out_file, n_threads, eval)
    return alignment

if __name__ == '__main__':
    main()
