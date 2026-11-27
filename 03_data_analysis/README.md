# Data Analysis

Example of the scripts used to process and analyze the simulation data.
Excerpts of the simulation outputs are given as examples to run these scripts.

Please note that averaging the combination of these scripts and data were not meant to be able to reproduce the same results as shown in the publication.
These are simplified Python notebooks featuring limited data, generally on only one replicate.
They still were very closely based on the scripts I used to produce the results shown in the publication.

Some explanation are also given on how they work and how the vast-pea objects and functions work.

## Files

### `scripts/`
- `0_preprocess/*`: Python scripts and notebooks on how to prepare the raw data into smaller and faster formats before proceeding with the analysis scripts
- `1_structure/*`: Python notebooks on how structure properties were calculated, including crosslink density, density, skewness-voxel, radial distribution function and conversion of curing simulations
- `2_thermomechanical/*`: Python notebooks on how thermomechanical properties were calculated, including Tg, coefficient of thermal expansion, free volume fraction and Young's modulus
- `3_dynamics/*`: Python notebooks on how dynamics were calculated, including number of BERs, bond autocorrelation function, mean square displacements, diffusion coefficients, segment orientation autocorrelation function and Arrhenian-like fits of bond and segmental lifetimes
- `fig_to_script.md`: index of which example scripts corresponds to which figure/table in the publication.

### `data/`
Excerpts of the simulation outputs.
Directories `0`, `1` and `2` represent my file structure for replicates. Only small portions of the thermo logs or dump trajectories are given.
Typically, only the data for the DGEBA/AFD is given, however, some properties required different processing for the different systems.
In these cases, data given and the processing is showcased accordingly.
- `BER/afd_ber_325K.bonds`: bond list dump file
- `BER/afd_ber_325K.dump`: trajectory dump file for the DGEBA/AFD system
- `BER/oaps_ber_325K.dump`: trajectory dump file for the DGEBA/AFD - OAPS system
- `BER/merge/*`: incomplete trajectory dump files used to demonstrate `scripts/0_preprocess/merge_trajectories.ipynb` where multiple dump files are merged into a single trajectory object
- `BER/save/afd_325K_all.traj_data.atoms.nc`: netCDF4 file used to restore a `trajectory` object containing all atoms.
    It is accompanied by a `<name>.details.dc` printed directory containing details about how the `trajectory` object was read and processed.
- `BER/save/afd_325K_t1-13.traj_data.atoms.nc`: netCDF4 file used to restore a `trajectory` object containing only atom types 1 and 13 (c and s1).
    Also has a `details.dc` file.
- `BER/save/afd_*K.bond_data.bonds.nc`: netCDF4 file used to restore a `bond_dynamics` object containing information about bond exchange reactions
- `BER/save/afd_*K.segment_data.*.nc`: netCDF4 files used to restore a `segment` object containing the trajectory of center of mass of segments, their directional vectors and information about beginnings and ends.
    Also has a `details.dc` file
- `dilatometry/output.tg`: LAMMPS thermo log output for the dilatometry protocol
- `dilatometry/system_tg.600.data`: LAMMPS data file obtained from the dilatometry protocol.
    Used to calculate FFV
- `dilatometry/ffv_logs/*.log.lunar`: LUNAR outputs from the calculation of FFV
- `young/afd_young.log.lammps`: LAMMPS thermo log output for the RFR straining protocol
- `preparation/afd_annealing2.data`: LAMMPS data file obtained after the second annealing protocol for the DGEBA/AFD system.
    Used to calculate crosslink density
- `preparation/oaps_annealing2.data`: LAMMPS data file obtained after the second annealing protocol for the DGEBA/AFD - OAPS system.
    Used to calculate crosslink density
- `preparation/*/afd_annealing.output`: LAMMPS thermo log output of the annealing equilibration for each replicate
- `preparation/*/afd_curing.output`: LAMMPS thermo log output of the curing simulation for each replicate

## Execution
No order is needed for the execution of the example scripts.
All raw or saved files are given to allow straightforward execution.
However, there is an order typically used for producing data from scratch:

### For bond dynamics (in order):
- `0_preprocess/extract_bonds.py`: Read bond list and extract bond dynamics object to saved `.bond_data.bonds.nc` file
- Compute properties (in no specific orders):
    - `3_dynamics/number_bers.ipynb`: Number of BERs
    - `3_dynamics/bond_autocorrelation.ipynb`: Bond autocorrelation function
    - `3_dynamics/arrhenian_fits.ipynb`: Arrhenian-like fit of bond lifetimes

### For trajectories (in order):
- `0_preprocess/extract_trajectory.py`: Read LAMMPS dump trajectory file and extract it to saved `.traj_data_atoms.nc`.
    Typically both all atoms and filtered for types c and s1 are extracted
- Compute properties (in no specific order):
    - `1_structure/rdf.ipynb`: Radial distribution function
    - `1_structure/skewness_voxel.ipynb`: Skewness-voxel
    - `3_dynamics/msd.ipynb`: Mean square displacements
    - `3_dynamics/diffusion.ipynb`: Diffusion coefficients

### For segments (in order):
- `0_preprocess/extract_trajectory.py`: Read LAMMPS dump trajectory file and extract it to saved `.traj_data_atoms.nc`.
    Extract for all atoms
- `0_preprocess/extract_segments.py`: Read the saved `trajectory` object, calculate segments and extract them to `.segment_data.coms.nc`, `.segment_data.ends.nc` and `.segment_data.vectors.nc` files
- Compute properties (in no specific order):
    - `3_dynamics/segment_orientation_autocorrelation.ipynb`: Segment orientation Legendre polynomial
    - `3_dynamics/arrhenian_fits.ipynb`: Arrhenian-life fit of segment lifetimes

### For Fractional Free Volume (in order):
- `0_preprocess/extract_ffv.ipynb`: Read LAMMPS data file and calculate the fractional free volume using LUNAR
- `2_thermomechanical/ffv.ipynb`: Read the LUNAR log output to obtain the fractional free volume

## Requirements
- Python 3.10+: Began this work on 3.10.8, currently on 3.12.3.
- [LUNAR](https://github.com/CMMRLab/LUNAR): version used is 21July2025
- [RFR](https://github.com/CMMRLab/RFR_stress_strain/tree/main): version used is 24March2025
- All non-default Python packages and the version used
    - lammps-logfile: 1.1.3
    - matplotlib: 3.11.2
    - numpy: 2.4.2
    - pandas: 3.0.1
    - Pint: 0.25.2
    - scipy: 1.17.0
    - xarray: 2026.2.0

## Todo
- [ ] Conversion doesn't feature convergence fit