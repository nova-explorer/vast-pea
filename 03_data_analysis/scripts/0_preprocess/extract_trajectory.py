#!/bin/python3

"""
Script to extract the data from trajectory dumps and save them to xarray netcdf files.
Those save files are smaller and faster to read, being ideal for large dumps.

This would typically be run in a loop over temperatures and replicates.

See 1_structure or 3_dynamics scripts for how to restore and use those saved files

Requirements:
    - lammps_traj from the vast_pea package
    - xarray
    - re
    - pandas
    - numpy
    - glob

Run using:
    ./extract_trajectory.py

Outputs save files
"""
from vast_pea.lammps_traj import trajectory

#########################################
# Reading and saving select atoms types #
#########################################

# Read the trajectory
traj = trajectory(path="../../data/BER/", # path to file
                  pattern="afd_ber_325K.dump", # LAMMPS trajectory file
                  filter_type="include", # "include" selects only atom types matching those in filter and "exclude" selects only atom types not matching those in filter. False can be used to disable the filter and keep all atoms as demonstrated later
                  filter=[1,13], # atom types to keep or remove. If filter_type=False, then you need to filter=False as well. In this case, it corresponds to the center carbons of DGEBA and sulfur atoms.
                  restore=False, # whether or not the trajectory is read from a saved netCDF file
                  verbose=True # print information
                  )

# Save the trajectory to a netCDF file
traj.save(name="afd_325K", # base name of the saved file. Filters and extension (.traj_data.atoms.nc) will be added automatically
          path="../../data/BER/save/", # path to save the file
          override=False # can be used to bypass the filter naming. Just change False to a string
          )

traj=None # for larger files, it can be useful to clear memory


######################################
# Reading and saving all atoms types #
######################################
# By setting filter_type and filter as False, all atoms are saved to the netCDF file.
# This however can take a lot of time and memory for long trajectories

traj = trajectory(path="../../data/BER/",
                  pattern="afd_ber_325K.dump",
                  filter_type=False,
                  filter=False,
                  restore=False,
                  verbose=True
                  )

# Since we don't have a filter, 'all' will be used instead of the filter types.
traj.save(name="afd_325K",
          path="../../data/BER/save/",
          override=False
          )