#!/bin/python3

"""
Script to extract the data from trajectory dumps and save them to xarray netcdf files.
Those save files are smaller and faster to read, being ideal for large dumps.

This would typically be run in a loop over temperatures and replicates.

Requirements:
    - lammps_traj found in src
    - xarray
    - re
    - pandas
    - numpy
    - glob

Run using:
    ./extract_trajectory.py

Outputs save files
"""
