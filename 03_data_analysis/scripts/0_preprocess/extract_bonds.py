#!/bin/python3

"""
Script to extract the data from bond dumps and save them to xarray netcdf files.
Those save files are smaller and faster to read, being ideal for large dumps.

This would typically be run in a loop over temperatures and replicates.

See 3_dynamics scripts for how to restore and use those saved files

Requirements:
    - lammps_bonds from the vast_pea package
    - xarray
    - pandas
    - numpy
    - glob

Run using:
    ./extract_trajectory.py

Outputs save files
"""
from vast_pea.lammps_bonds import bond_dynamics

#########################################
# Reading and saving select atoms types #
#########################################

# Read the bonds
bd = bond_dynamics(path="../../data/BER/", # path to file
                   pattern="afd_ber_325K.bonds", # LAMMPS trajectory file
                   restore=False, # whether or not the trajectory is read from a saved netCDF file
                   verbose=True # print information
                   )

# Save the bonds to a netCDF file
bd.save(name="afd_325K", # base name of the saved file. Filters and extension (.traj_data.atoms.nc) will be added automatically
        path="../../data/BER/save/", # path to save the file
        )
