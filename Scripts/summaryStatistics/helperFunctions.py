"""
This module contains general functions for working with summary statistic data
Import modules and add functions as necessary.
"""

import os
import pickle
import re #read regex
    
def create_pickle_list(directory):  
    """
    Creates a list of pickle file names.
    
    @param  directory   string containing .pickle files
    @return pickle_list list of .pickle files
    """
    pickle_list = []
    for file in sorted(os.listdir(directory)):
        if file.endswith(".pickle"):
            pickle_list.append(file)
    return pickle_list
