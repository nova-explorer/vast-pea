# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: wrapper for the lammps_logfile package

TODO:
    -
"""

import pandas as pd
import numpy as np
import lammps_logfile as ll
import os
from typing import List, Tuple

def is_line_all_number(line: str) -> bool:

    is_first_number = line.strip()[0] in '0123456789'
    is_no_char = not any(elem in line for elem in 'abcdfghijklmnopqrstuvwxyz') # removed 'e' as it is part of scientific notation used by LAMMPS

    if is_first_number and is_no_char:

        return True

    else:

        return False

def read_lines(filename: str) -> List[str]:

    with open(filename, 'rt') as f:
        lines = f.readlines()

    return lines

def preprocess_lines(lines: List[str]) -> List[str]:
    total_lines = len(lines)

    lastkeywords = ""
    lastfeeter = ""
    data = []
    feeter = ""

    for i in range(total_lines):

        if "Step" in lines[i]:

            header = lines[i-1]
            keywords = lines[i].strip()

            feeter, counter = find_next_feeter(lines, i, total_lines)

            if keywords != lastkeywords or counter > 5:

                data.append(lastfeeter)
                data.append(header)
                data.append(keywords + "\n")

                lastkeywords = keywords
                lastfeeter = feeter

            data += [ "\t" + k.strip() + "\n" for k in lines[i+1:i+counter+1] if "WARNING" not in k ]

    data.append(feeter)

    return data

def find_next_feeter(lines: List[str], i: int, total: int) -> Tuple[str, int]:

    counter=0
    feeter = ""

    for j in range(i+1,total):

        if "WARNING" in lines[j]:
            continue

        else:

            if lines[j].strip() and is_line_all_number(lines[j]):
                counter+=1

            else:
                feeter = lines[j]
                break

    return feeter, counter

def write_tmpfile(lines: List[str]) -> None:
    with open(os.getcwd()+"/tmpfile", 'wt') as f:
        f.writelines(lines)

def remove_tmpfile() -> None:
    os.remove(os.getcwd()+"/tmpfile")

def thermo_to_df(thermo: ll.File, skip_first: int, skip_last: int, timestep: float) -> pd.DataFrame:
    cols = thermo.get_keywords()
    nb_runs = thermo.get_num_partial_logs()
    d = {}

    for col in cols: # pyright: ignore[reportOptionalIterable]
        d.update({col: special_get(thermo, col, nb_runs, skip_first, skip_last)})
    df = pd.DataFrame(data=d, index=d['Step'])

    df = df.drop(["Step"], axis=1)

    df.index.name = "Step"
    df.index.astype(int)

    df.index = df.index / 1e6 * timestep

    return df

def special_get(thermo: ll.File, col: str, nb_runs: int, skip_first: int, skip_last: int) -> np.ndarray:
    array = np.array([])

    for i in range(nb_runs+1)[:-1]:
        data = thermo.get(col, run_num=i)

        if data is None:
            array = np.concatenate(( array, [None]*len(thermo.get("Step", run_num=i)[skip_first:-skip_last or None]) )) # type: ignore

        else:
            array = np.concatenate(( array, data[skip_first:-skip_last or None] ))

    return array

def lammps_out_to_pd(path: str, pattern: str, n_replicas: int, skip_first: int=50, skip_last: int=50, timestep: float=1, timestep_offset=[0]) -> pd.DataFrame:

    dfs = pd.DataFrame()

    if len(timestep_offset) == 1:
        timestep_offset = timestep_offset * n_replicas

    for i in range(n_replicas):

        if n_replicas > 1:
            filename = "{}/{}/{}".format(path, str(i), pattern)

        else:
            filename = "{}/{}".format(path, pattern)

        lines = read_lines(filename)
        lines = preprocess_lines(lines)

        write_tmpfile(lines)
        thermo = ll.File(os.getcwd()+"/tmpfile")
        remove_tmpfile()

        df = thermo_to_df(thermo, skip_first, skip_last, timestep)

        df.index += timestep_offset[i]

        df = df.sort_index()

        df['replica'] = int(i)

        dfs = pd.concat([dfs, df])

    return dfs