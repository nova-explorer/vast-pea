# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: trajectory class used to read and process LAMMPS dump trajectory files

TODO:
    - select n frames
    - set columns type properly
"""

import xarray as xr
import re
import pandas as pd
import numpy as np

import io_

class trajectory():
    def __init__(self,
                 path = "./", pattern = "*.dump",
                 filter_type = "include", filter = False,
                 restore = False, verbose = True
    ) -> None:
        self.verbose = verbose

        self.file_list = io_.create_file_list(path, pattern)

        self.__print( "\tFile list contains {} files\n".format(len(self.file_list)) )

        self.filter_type = filter_type
        self.filter = filter

        if not restore:

            self.__print("\tReading trajectory...\n")
            self._atoms = self.__read_dumps()

        else:

            self.__restore()

        self.nb_atoms = len(self._atoms.id.values)
        self.timesteps = self._atoms.ts.values
        self.ids = self._atoms.id.values
        self.comps = self._atoms.comp.values

        self.__print("\tTrajectory done!\n")

    def __read_dumps(self) -> xr.Dataset:
        comp = ['x', 'y', 'z']
        nb_atoms = 0
        step = 0
        bounds_array = []
        timestep = []
        data = []

        for filename in self.file_list:
            lines = io_.read_file(filename)
            cols, header = self.__get_dump_cols(lines)

            for i, line in enumerate(lines):

                if "ITEM: TIMESTEP" in line:
                    step = int(lines[i+1])

                elif "ITEM: NUMBER OF ATOMS" in line:
                    nb_atoms = int(lines[i+1])

                elif "ITEM: BOX BOUNDS pp pp pp" in line:
                    bounds = []
                    for j in range(len(comp)):
                        bounds.append([ float(val) for val in lines[i+j+1].strip().split() ])
                    bounds_array.append([ j[1] - j[0] for j in bounds ])

                elif header in line:
                    data.append(self.__lines_to_df(lines[i+1:i+nb_atoms+1], cols))
                    timestep.append(step)

        data = self.__dfs_to_ds(data, timestep)
        data['bounds'] = xr.DataArray(bounds_array, coords = [data.ts, comp], dims = ['ts', 'comp'])

        return data

    def __get_dump_cols(self, lines:list) -> list:
        counter = 0
        while not re.search("ITEM: ATOMS", str(lines[counter])):
            counter += 1
        return lines[counter].strip().split() [2:], lines[counter].strip() # pyright: ignore[reportReturnType]

    def __lines_to_df(self, lines, column_names) -> pd.DataFrame:
        data = []
        type_id = None
        # Checks the position of the property: type.
        for i, name in enumerate(column_names):
            if name == 'type':
                type_id = i

        if 'xu' in column_names:
            column_names[column_names.index('xu')] = 'x'
        if 'yu' in column_names:
            column_names[column_names.index('yu')] = 'y'
        if 'zu' in column_names:
            column_names[column_names.index('zu')] = 'z'

        match self.filter_type:
            case False:
                for i in lines:
                    data.append(i.strip().split())

            case "exclude":
                for i in lines:
                    line = i.strip().split()
                    if int(line[type_id]) not in self.filter: # pyright: ignore[reportOperatorIssue]
                        data.append(i.strip().split())

            case "include":
                for i in lines:
                    line = i.strip().split()
                    if int(line[type_id]) in self.filter: # pyright: ignore[reportOperatorIssue]
                        data.append(i.strip().split())

        df = pd.DataFrame(data, columns = column_names)
        df = df.apply(pd.to_numeric)
        df = df.set_index(['id'])
        df = df.sort_index()

        return df

    def __dfs_to_ds(self, dfs, timesteps) -> xr.Dataset:
        ds = []
        for i, df in enumerate(dfs):
            ds.append(xr.Dataset())
            for j in df:
                ds[i][j] = df[j].to_xarray()
            ds[i].update( {'ts': ( 'ts', [timesteps[i]] )} )
        return xr.concat(ds, dim = 'ts', data_vars='all')

    def __print(self, text) -> None:
        if self.verbose:
            print(text, end="\r")

    def __restore(self) -> None:
        CURRENT_MAX_FILES = 1
        currently_imported = {'atoms':False}

        if len(self.file_list) > CURRENT_MAX_FILES:
            raise EnvironmentError("Too many files for trajectory restore: "+", ".join(self.file_list))

        for filename in self.file_list:
            if ".traj_data.atoms.nc" in filename and not currently_imported['atoms']:
                self._atoms = io_.read_dataset(filename)
                currently_imported['atoms'] = filename
            else:
                raise EnvironmentError("Already imported atoms: " + currently_imported['atoms'] + ", current file: " + filename) # pyright: ignore[reportOperatorIssue]

    def save(self, name, path, override=False) -> None:
        if override:
            filename_out = "{}_{}".format(name,
                                            override # pyright: ignore[reportGeneralTypeIssues]
                                            )
        elif self.filter:
            filename_out = "{}_t{}".format(name,
                                            "-".join([str(i) for i in self.filter]) # pyright: ignore[reportGeneralTypeIssues]
                                            )
        else:
            filename_out = "{}_all".format(name)

        details = {'file_list': self.file_list,
                   'filter_type': self.filter_type,
                   'filter': self.filter
                   }

        io_.save_xarray(self._atoms,
                        path,
                        filename_out+".traj_data.atoms.nc"
                        )

        io_.save_dict(details,
                      path,
                      filename_out+".traj_data.details.dc"
                      )

    def get_atoms(self) -> xr.Dataset:
        return self._atoms

    def set_atoms(self, new_atoms) -> None:
        self._atoms = new_atoms.copy()
        self.nb_atoms = len(self._atoms.id.values)
        self.ids = self._atoms.id.values
        self.timesteps = self._atoms.ts.values

    def get_unwrapped(self) -> xr.DataArray:
        return self.get_positions() + self.get_images() * self.get_bounds()

    def get_positions(self) -> xr.DataArray:
        return xr.DataArray(np.transpose([self._atoms.x, self._atoms.y, self._atoms.z]), coords=[self.ids, self.timesteps, self.comps], dims=['id', 'ts', 'comp'])

    def get_images(self) -> xr.DataArray:
        return xr.DataArray(np.transpose([self._atoms.ix, self._atoms.iy, self._atoms.iz]), coords=[self.ids, self.timesteps, self.comps], dims=['id', 'ts', 'comp'])

    def get_bounds(self) -> xr.DataArray:
        return self._atoms.bounds

    def add_bond_vectors(self, bond_connections) -> xr.Dataset:

        total_ts = len(bond_connections)
        total_ids = len(bond_connections[0])

        if total_ts != len(self.timesteps):
            raise Warning("Different number of timesteps between atoms and bonds trajectories")
        # if total_ids*2 != self.nb_atoms:
        #     raise Warning("Different number of particles between atoms and bonds trajectories")

        coords = ['x', 'y', 'z']
        icoords = ['ix', 'iy', 'iz']
        # creates a list containing all the columns that are not positions or images
        droppers = []
        for i in self._atoms.data_vars:
            if i not in coords and i not in icoords:
                droppers.append(i)
        data = self._atoms.drop_vars(droppers)

        vectors_traj = []
        positions_traj = []
        images_traj = []

        #probably remove those once I'm done testing stuff
        ids_traj = []
        pos1_traj = []
        pos2_traj = []

        for ts in range(total_ts):

            atom0 = data.sel(ts=self.timesteps[ts], id=list(bond_connections[ts].keys()))
            atom0_images = np.array([atom0.ix.values, atom0.iy.values, atom0.iz.values])
            atom0 = np.array([atom0.x.values, atom0.y.values, atom0.z.values])

            atom1 = data.sel(ts=self.timesteps[ts], id=list(bond_connections[ts].values()))
            atom1 = np.array([atom1.x.values, atom1.y.values, atom1.z.values])

            direction = np.where(atom0 < atom1, -1, 1)
            bounds = self._atoms.bounds.isel(ts=ts).values

            for i in range(len(coords)):
                atom1[i] = atom1[i] + direction[i] * bounds[i] * np.round((atom0[i]-atom1[i])/bounds[i])

            v = atom1 - atom0
            v /= np.sqrt((v**2).sum(axis=0))

            cm = (atom0 + atom1)/2

            vectors_traj.append(v)
            positions_traj.append(cm)
            images_traj.append(atom0_images)
            ids_traj.append(np.array([
                list(bond_connections[ts].keys()),
                list(bond_connections[ts].values())
                ]).transpose()
            )
            pos1_traj.append(atom0)
            pos2_traj.append(atom1)

        ds = xr.Dataset(
            data_vars={
                'vector':(['ts','comp','id'], vectors_traj),
                'positions':(['ts','comp','id'], positions_traj),
                'images':(['ts','comp','id'], images_traj),
                'atom_ids':(['ts', 'id', 'atom'], ids_traj),
                'pos1':(['ts','comp', 'id'], pos1_traj),
                'pos2':(['ts', 'comp', 'id'], pos2_traj)
                },
            coords={'ts':self.timesteps, 'comp':coords, 'id':range(total_ids), 'atom':['0','1']}
        )

        return ds.transpose('ts', 'id', 'comp', 'atom')

    def adjust_molids(self, new_molids: xr.DataArray) -> None:
        # new_molids should only be function of ids, not comp/timestep
        assert len(new_molids.shape) == 1

        # check if atoms and new_molids have same ids
        assert sorted(self.ids) == sorted(new_molids.id.values)

        df = self._atoms.to_dataframe()

        for i in self.ids:
            df.loc[(slice(None),i),'mol'] = new_molids.sel(id=i)

        self._atoms = df.to_xarray()

        self._atoms['mol'] = self._atoms['mol'].mean(dim='comp').astype(int)
        self._atoms['type'] = self._atoms['type'].mean(dim='comp').astype(int)
        self._atoms['x'] = self._atoms['x'].mean(dim='comp')
        self._atoms['y'] = self._atoms['y'].mean(dim='comp')
        self._atoms['z'] = self._atoms['z'].mean(dim='comp')
        self._atoms['ix'] = self._atoms['ix'].mean(dim='comp').astype(int)
        self._atoms['iy'] = self._atoms['iy'].mean(dim='comp').astype(int)
        self._atoms['iz'] = self._atoms['iz'].mean(dim='comp').astype(int)
        self._atoms['bounds'] = self._atoms['bounds'].mean(dim='id')

    def agglomerate(self) -> None:

        positions = self.get_positions().to_numpy()

        for comp_id in range(len(self.comps)):
            comp = self.comps[comp_id]

            dist = positions[:,:,comp_id] - [positions[0,:,comp_id]]

            bound = self._atoms.bounds.sel(comp=comp).to_numpy()

            self._atoms[comp] -= (bound * np.round(dist/bound)).T

    def recenter(self, masses) -> None:
        mass_sum = masses.sum().values
        for i in self.comps:
            com = (self._atoms[i] * masses).sum(dim='id') / mass_sum
            self._atoms[i] -= com

    def reset_images(self) -> None:
        for i in self.comps:
            self._atoms['i'+i] *= 0
