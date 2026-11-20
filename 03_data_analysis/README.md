# Data Analysis

Some details

## Files
- `src/*.py`: backend classes and functions used by the `scripts/` frontend files
- `data/*`: excerpts of the simulation outputs. Directories `0`, `1` and `2` represent my file structure for replicates. Only small portions of the thermo logs or dump trajectories are given. Typically, only the data for the DGEBA/AFD is given, however, some properties required different processing for the different systems. In these cases, data given and the processing is showcased accordingly.
- list of files

## Execution
- create saves

- what each file does

## Requirements
- Python 3.10+: Began this work on 3.10.8, currently on 3.12.3.
- Packages

## TODO
- [ ] comment src/
    - [ ] computers.py
    - [x] constants.py
    - [ ] fitting.py
    - [ ] io_.py
    - [ ] lammps_bonds.py
    - [ ] lammps_data.py
    - [ ] lammps_out.py
    - [ ] lammps_segments.py
    - [x] masses.py
    - [x] processing.py
- [ ] package src instead of adding to path
- [ ] consider changing crosslink_density.ipynb trajectory read to restore
- [ ] consider adding merge_trajectories for run with multiple dump files