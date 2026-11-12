#!/bin/python3

"""
Recenters the monomer agglomeration and removes excess volume

Needs the data_file class which can be found in 03_data_analysis/src/lammps_data

Run using:
    ./recenter.py system_equil-small1.data

Outputs system_equil-small1_centered.data
"""
import sys
import os

from lammps_data import data_file

def main(input_file, path):
    data = data_file(path=path, pattern=input_file)

    data.agglomerate()
    data.recenter()
    data.reset_images()

    output_file = path + "".join(input_file.split('.')[:-1]) + "_centered.data"

    data.write_file(output_file)

if __name__ == "__main__":

    input_file = sys.argv[1]
    path = os.getcwd() + "/"

main(input_file, path)