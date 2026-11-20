# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: Mass for atoms types in PCFF that I wanted to have quick reference
"""

# Mass of atom types in PCFF
PCFF_MASSES = {
    'c'  : 12.01115,
    'c1' : 12.01115,
    'c2' : 12.01115,
    'c3' : 12.01115,
    'cp' : 12.01115,
    'h'  : 1.00797,
    'hc' : 1.00797,
    'hn' : 1.00797,
    'ho2': 1.00797,
    'nn' : 14.0067,
    'o=' : 15.9994,
    'oc' : 15.9994,
    'oh' : 15.9994,
    'osi': 15.9994,
    's1' : 32.0640,
    'sf' : 32.0640,
    'sio': 28.0860,
}

def types_to_masses(type_list: list[str], mass_dict: dict[float]) -> list[float]:
    """Converts a list of atom types to their mass equivalent.

    Args:
        type_list (list[str]): list of atoms types in the system ordered by their LAMMPS id.
        mass_dict (dict[float]): PCFF_MASSES like dict.

    Returns:
        list[float]: list of masses of the atoms types in the system ordered by their LAMMPS id.
    """

    return [mass_dict[i] for i in type_list]