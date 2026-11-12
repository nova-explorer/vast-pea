# System Creation

This is where the initial systems were generated.

The molecular structures were created using [Avogadro2](https://avogadro.cc/index.html) and are in the `0_molecules` directory.

The [LUNAR](https://github.com/CMMRLab/LUNAR) toolkit was then used to transfer those to LAMMPS data files, merge forcefield parameters with the pre/post-reaction templates and pack the initial systems which were then used for simulations.

## Files
- `0_molecules/*.mol2` files: Structures of the molecules and pre/post-reaction templates created in Avogadro2.

- `creator.ipynb`: Wrapper script for LUNAR. It was made mostly by taking the parts of the LUNAR scripts I needed and putting them together so I can Run All cells easily.
    Calculates automatically the duplicate variable in cell_builder so that the system created matches as closely as possible a desired number of atoms.
    Also runs cell_builder multiple times to create different initial configurations at once.

    It is located in the main LUNAR directory during execution.

- `PCFF.py`: Minor fix of the `LUNAR/src/atom_typing/PCFF.py` file to attribute types correctly for sio in POSS and lone hydrogens in post-reaction templates.

## Execution
- Move `creator.ipynb` to the main LUNAR directory
- Move `PCFF.py` to `LUNAR/src/atom_typing/`
- Run cells in creator.ipynb
- When asked for a directory, select the parent one to `0_molecules/`
- LUNAR will execute and create its typical directories, the built systems will be in `4_cell_builder/`, with each one in its separate directory

## Requirements
- LUNAR version 21July2025: Could probably run higher or lower version, but this is the one that was used
- Python 3.9+: Version used for this example was 3.10.12, but 3.12.3 worked as well
    - tkinter: only used for filedialog, could be replaced by manually specifying the parent directory of `0_molecules/`
    - time: only used to provide different seeds for cell_builder