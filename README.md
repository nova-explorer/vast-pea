# Vitrimer Atomistic Simulation Tools for Prediction and Exploration Analysis (VAST-PEA)

This is a collection of scripts, tools utilities and examples for the atomistic simulation of vitrimers.

These resources were written by Olivier Couture in the context of his PhD at the Luxembourg Institute of Science and Technology, and used in the publication **ADD LINK**.

The framework shown in that publication (and here) was used to generate initial systems (using [LUNAR](https://github.com/CMMRLab/LUNAR)), run the [LAMMPS](https://www.lammps.org) simulations, and process the data to study the structure, thermomechanical properties and dynamics of the simulations.

- [`01_system_creation/`](01_system_creation/README.md) details how systems were initially created using LUNAR
- [`02_simulations/`](02_simulations/README.md) lists the different LAMMPS simulations that were realized and in which order
- [`03_data_analysis/`](03_data_analysis/README.md) gives scripts and data examples calculating the different results shown in the publication.
    Do note that this directory features 92 files at 175MB, being much heavier that the rest of this repo
- [`vast_pea/`](vast_pea/) contains the code for the `vast_pea` Python package, acting essentially as a backend for the scripts in `03_data_analysis/`

The `vast_pea` package should be able to install with a simple
```
pip install .
```

This collection of LAMMPS and Python scripts is given openly in order to promote open science and potentially help others on their research journey.
They are, however, flawed and not completely documented (namely the `vast_pea` package).
They are also provided without guarantee of them working in your case, or promise that they will be maintained over time.

## Requirements
These are less requirements than what was used to develop and run the resources provided.

### System Creation
- LUNAR (version 21July2025)
- Python (3.12.3)

### Simulations
- LAMMPS (2Aug2023_update2-foss-2023a-kokko)
- Python (3.6.8)
- vast_pea package

### Data analysis
- LUNAR (version 21July2025)
- vast_pea package

### vast_pea package
- Python (3.12.3)
- lammps-logfile (1.1.3)
- matplotlib (3.11.2)
- numpy (2.4.2)
- pandas (3.0.1)
- Pint (0.25.2)
- scipy (1.17.0)
- xarray (2026.2.0)

## Terms of use
- Read and respect the [license](LICENSE)
- Acknowledge the use of any reports or publications of results obtained using the resources provided here by citing **ADD LINK**

## Acknowledgements
This collection (and the aforementioned publication) uses several toolkits/codes/modules that are worth acknowledging, since without their open source nature and availability, this project would be far behind what it is now. For this I thank the authors of :
- [LAMMPS](https://www.lammps.org): Simulation code
- [REACTER](https://www.reacter.org/): LAMMPS module enabling heuristic topology changes (reactions).
    Their [examples](https://github.com/lammps/lammps/tree/develop/examples/PACKAGES/reaction) were really useful in understanding how to properly structure a LAMMPS script.
- [LUNAR](https://github.com/CMMRLab/LUNAR): Toolkit for generating REACTER files and computing certain properties (free volume in this case).
- [RFR-method](https://github.com/CMMRLab/RFR_stress_strain): A useful method for easily analyzing stress-strain curves.
    Their LAMMPS scripts and Python code were adapted to compute the Young's modulus using their method.

## TODO
- [x] package src instead of adding to path

### simulations
- [ ] comment charge_rescale.py

### data_analysis
- [ ] consider changing crosslink_density.ipynb trajectory read to restore
- [x] consider adding merge_trajectories for run with multiple dump files

### vast-pea Python package
- [ ] comment vast_pea/
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
- [ ] Some instructions about usage and installation
- [ ] Add a better list of dependencies
- [ ] pattern should be changed to filename for functions that can only read 1 file
