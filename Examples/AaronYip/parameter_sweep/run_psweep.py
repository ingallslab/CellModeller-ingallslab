#!/usr/bin/env python3

import re
import os
import subprocess
import copy
import psweep as ps

# CellModeller module name
module_name = 'tutorial1.py'

def func(pset):
    """
    This function runs the simulation and is looped in ps.run_local()
    @param  pset    dict containing parameters and values (e.g., {"x": 1, "y": 2}) 
    @return cmd     the commands ran in the terminal   
    """
    export_path = os.path.join(pset["_calc_dir"], pset["_pset_id"])
    parameters = get_params_dict(pset)
    search_and_replace_parameters(module_name, parameters)
    cmd = (
        f"mkdir -p $(dirname {export_path}); "
        f"python batch_psweep.py {module_name} {export_path}"
    )
    subprocess.run(cmd, shell=True)
    return {"cmd": cmd}
    
def search_and_replace_parameters(file_path, parameters):
    """
    Search and replace parameter values in a file with the pattern:
    
    name = value
    
    @param  file_path   str of file name containing parameters to be modified
    @param  parameters  dict containing parameter names and values {name: value}
    """

    with open(file_path, "r+") as file:
        file_contents = file.read()
        for name,value in parameters.items():
            search_string = f'{name} = (\d+)'
            substitute_string = f'{name} = {value}'
            file_contents = re.sub(search_string, substitute_string, file_contents)
        file.seek(0)
        file.truncate()
        file.write(file_contents)


def get_params_dict(pset):
    """
    Creates a copy of pset with only the parameter values.
    A helper function for search_and_replace_parameters.
    @param  pset        dict containing parameters and values (e.g., {"x": 1, "y": 2})
    @return new_dict    pset dict without keys that are used for parsing through parameter sweep database
    """
    new_dict = copy.deepcopy(pset) #creates a separate copy of pset dect
    ignore_keys = ['_run_id', '_pset_id', '_calc_dir', '_time_utc', '_pset_sha1', '_pset_seq', '_run_seq']
    for key in ignore_keys:
        del new_dict[key]
    return new_dict


if __name__ == "__main__":
    params = ps.plist("growthRate", [1,2,3])
    df = ps.run_local(func, params)
    #print(df.columns)
