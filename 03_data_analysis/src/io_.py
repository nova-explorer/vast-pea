# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: Functions to read/write files

TODO:
    -
"""

import os
from glob import glob
import xarray as xr
import numpy as np

def create_file_list(path:str, file_pattern:str) -> list:
    file_list = glob(os.path.join(path, file_pattern))
    if len(file_list) == 0:
        raise EnvironmentError("No file matching found :"+os.path.join(path, file_pattern))
    return file_list

def make_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print("Could not create dir:", path)

def read_file(filename:str) -> list:
    with open(filename, 'rt') as f:
        data = f.readlines()
    return data

def save_xarray(ds, path, name, format_type=np.float32):
    """
    TODO: convert to os.joinpath if path issues
    """
    make_dir(path)
    ds = ds.astype(dtype=format_type)
    ds.to_netcdf(path+name)

def save_dict(dict, path, name):
    make_dir(path)
    with open (path+name, 'wt') as f:
        f.write(str(dict))

def read_dataset(filename):
    return xr.open_dataset(filename)

def read_dataarray(filname):
    return xr.open_dataarray(filname)

def to_dump(data, bounds, name, path):
    filename = os.path.join(path, name)
    line_template = "{} 1 {} {} {}\n"

    with open(filename, 'wt') as f:
        for ts in data.ts.values:
            f.write( 'ITEM: TIMESTEP\n' )
            f.write( str(ts)+'\n' )

            f.write( 'ITEM: NUMBER OF ATOMS\n' )
            f.write( str(len(data.id.values)) + '\n' )

            f.write( 'ITEM: BOX BOUNDS pp pp pp\n' )
            for comp in bounds.comp.values:
                f.write('0 ' + str(bounds.sel(ts=ts, comp=comp).values) + '\n')

            f.write( 'ITEM: ATOMS id type x y z\n')
            for i in data.id.values:
                f.write(line_template.format(i,
                                             *data.sel(ts=ts, id=i).values
                                             )
                                             )

def combine_to_dump(atoms, cms, name, path):
    filename = os.path.join(path, name)
    atoms_template = "{id} 1 {type} {x} {y} {z} {ix} {iy} {iz}\n"
    cms_template = "{id} 1 {type} {x} {y} {z} 0 0 0\n"
    cms_type = atoms.type.max().values+1

    n_all = len(atoms.id.values) + len(cms.id.values)
    n_atoms = len(atoms.id.values)
    with open(filename, 'wt') as f:
        for ts in atoms.ts.values:
            f.write( 'ITEM: TIMESTEP\n' )
            f.write( str(ts)+'\n' )

            f.write( 'ITEM: NUMBER OF ATOMS\n' )
            f.write( str(n_all) + '\n' )

            f.write( 'ITEM: BOX BOUNDS pp pp pp\n' )
            for comp in atoms.bounds.comp.values:
                f.write('0 ' + str(atoms.bounds.sel(ts=ts, comp=comp).values) + '\n')

            f.write( 'ITEM: ATOMS id mol type x y z ix iy iz\n')
            for i in atoms.id.values:
                f.write(atoms_template.format(
                    id=i,
                    type=atoms.sel(ts=ts, id=i).type.values,
                    x=atoms.sel(ts=ts, id=i).x.values,
                    y=atoms.sel(ts=ts, id=i).y.values,
                    z=atoms.sel(ts=ts, id=i).z.values,
                    ix=atoms.sel(ts=ts, id=i).ix.values,
                    iy=atoms.sel(ts=ts, id=i).iy.values,
                    iz=atoms.sel(ts=ts, id=i).iz.values,
                    )
                    )

            for i, id_i in enumerate(cms.id.values):
                f.write(cms_template.format(
                    id=i + n_atoms + 1,
                    type=cms_type,
                    x=cms.sel(ts=ts, id=id_i, comp='x').values,
                    y=cms.sel(ts=ts, id=id_i, comp='y').values,
                    z=cms.sel(ts=ts, id=id_i, comp='z').values
                    )
                    )

def to_data(filename, header, bounds, labels, masses, coeffs, atoms, velocities, interactions) -> None:
    bound_template = "{} {} {}lo {}hi\n"
    mass_template = "{}   {}  # {}\n"
    atoms_template = "  {id}     {mol}   {type}   {charge}   {x}    {y}    {z}   {ix}  {iy}  {iz}\n"

    with open(filename, 'wt') as f:

        # write header
        for line in header:
            f.write(line)

        # write bounds
        for comp in bounds:
                f.write(bound_template.format(
                                              bounds[comp][0],
                                              bounds[comp][1],
                                              comp,
                                              comp
                                              ))

        f.write('\n')

        # write labels
        for line in labels:
            f.write(line)

        # write masses
        f.write('Masses \n\n')

        for i, key in enumerate(masses):
            f.write(mass_template.format(i+1, masses[key], key))

        f.write('\n')

        # write coeffs
        for line in coeffs:
            f.write(line)

        # write atoms
        f.write('Atoms # full\n\n')

        for i in atoms.id.values:
            f.write(atoms_template.format(
                id=i,
                mol=   atoms.sel(id=i).molecule.values,
                type=  atoms.sel(id=i).type.values,
                charge=atoms.sel(id=i).charge.values,
                x=     atoms.sel(id=i).x.values,
                y=     atoms.sel(id=i).y.values,
                z=     atoms.sel(id=i).z.values,
                ix=    atoms.sel(id=i).ix.values,
                iy=    atoms.sel(id=i).iy.values,
                iz=    atoms.sel(id=i).iz.values
            ))

        f.write('\n')

        if velocities:
            for line in velocities:
                f.write(line)

        # write interactions
        for line in interactions:
            f.write(line)
