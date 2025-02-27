#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script contains general functions used in several scripts
within the pipeline. 

This file should not be run on its own as a script
"""
from os import listdir
from os.path import exists
import subprocess as sp

def file_path_getter(dir, extension = None):
    """This function get the paths to all files located
    in sub directories

    Args:
        dir (str): path to main directory which should contain sub directories
        extension (str): the extension of the file you want to select. Defaults to None.
    Returns:
        list: paths to all the files
    """
    files = []
    for sub_dir in listdir(dir):
        if sub_dir.find('.') == -1:
            for file in listdir(f'{dir}/{sub_dir}'):
                if extension != None and file.endswith(extension):
                    files.append(f'{dir}/{sub_dir}/{file}')
                elif extension == None:
                    files.append(f'{dir}/{sub_dir}/{file}')
    return files

def dir_maker(dir_path):
    """This function checks if the directory exists
    and if not it makes the directory

    Args:
        dir_path (str): path to where the directory should be written
    Returns:
        str: path to where the directory is written
    """
    if exists(dir_path) == False:
        sp.run(f'mkdir -p {dir_path}', shell = True)
    return dir_path