import xarray as xr
import pandas as pd
import numpy as np
import re

import io_

"""
TODO:
    - check if file splitting works correctly
    - and workshop options for very large files
    - better print
"""
class bond_dynamics():
    def __init__(self,
                 path = "./", pattern = "bonds.dump",
                 restore = False, verbose = True
                 ) -> None:
        self.verbose = verbose

        self.file_list = io_.create_file_list(path, pattern)
        self.__print( "\tFile list contains {} files\n".format(len(self.file_list)) )

        if not restore:

            self.__print("\tReading bond dynamics...\n")
            self._bonds = self.__read_bond_dumps()

        else:

            self.__restore()

        self.nb_bonds = len(self._bonds.id.values)
        self.timesteps = self._bonds.ts.values
        self.ids = self._bonds.id.values

        self.__print("\tBond Dynamics done!\n")

    def __read_bond_dumps(self) -> xr.Dataset:

        cols = ''
        nb = 0

        step = None
        timestep = []
        data = []

        for filename in self.file_list:
            lines = io_.read_file(filename)
            cols, header = self.__get_dump_cols(lines) # could run only one time

            for i, line in enumerate(lines):
                line = line.rstrip()

                if line == "ITEM: TIMESTEP":
                    step = int(lines[i+1])

                elif line == "ITEM: NUMBER OF ENTRIES":
                    nb = int(lines[i+1])

                elif header in line:
                    data.append(self.__str_to_list(lines[i+1:i+nb+1]))
                    timestep.append(step)

        data = self.__list_to_ds(data, timestep, cols)

        return data

    def __get_dump_cols(self, lines:list):
        counter = 0
        while not re.search("ITEM: ENTRIES", str(lines[counter])):
            counter += 1
        return lines[counter].strip().split() [2:], lines[counter].strip()

    def __str_to_list(self, lines) -> list:
        return [[int(j) for j in i.rstrip().split(' ')] for i in lines]

    def __list_to_ds(self, data, timestep, cols) -> xr.Dataset:
        data = np.array(data)
        bond_ids = range(len(data[0]))

        data_vars = { cols[i]: (['ts', 'id'], data[:,:,i]) for i in range(len(cols)) }

        ds = xr.Dataset(data_vars, coords={'ts':timestep, 'id':bond_ids})

        return ds

    def __print(self, text) -> None:
        if self.verbose:
            print(text, end="\r")

    def __restore(self) -> None:
        CURRENT_MAX_FILES = 1
        currently_imported = {'bonds':False}

        if len(self.file_list) > CURRENT_MAX_FILES:
            raise EnvironmentError("Too many files for bond restore: " + ", ".join(self.file_list))

        for filename in self.file_list:
            if ".bond_data.bonds.nc" in filename and not currently_imported['bonds']:
                self._bonds = io_.read_dataset(filename)
                currently_imported['bonds'] = filename
            else:
                raise EnvironmentError("Already imported bonds: " + currently_imported['bonds'] + ", current file: " + filename) # pyright: ignore[reportOperatorIssue]

    def save(self, name, path) -> None:
        filename_out = "{}".format(name)

        # not saving details since there's none relevent for now

        io_.save_xarray(self._bonds,
                        path,
                        filename_out+".bond_data.bonds.nc"
                        )

    def get_bonds_setlist(self) -> list:
        # TODO: optimize

        str_template = "{} {}"
        data = np.array(self.get_bond_atoms())

        data = [ set([ str_template.format(data[0,ts,bid], data[1,ts,bid]) for bid in range(self.nb_bonds) ]) for ts in range(len(self.timesteps)) ]

        return data

    def get_connections(self) -> list:

        connections = []

        for ts in range(len(self.timesteps)):
            data = self._bonds.isel(ts=ts)

            atom1 = data['c_batom1'].values.astype(int)
            atom2 = data['c_batom2'].values.astype(int)

            connections.append({atom1[i]:atom2[i] for i in range(self.nb_bonds)})

        return connections

    def get_bond_atoms(self) -> list:
        """
        Format of output \\
        [ \\
            [ \\
                [id1, id2], ts1 \\
                [...], ts2 \\
            ], atom1 \\
        ]
        """
        return [self._bonds['c_batom1'].values, self._bonds['c_batom2'].values]

    def get_bonds(self) -> xr.Dataset:
        return self._bonds

    def set_bonds(self, new_bonds) -> None:
        self._bonds = new_bonds.copy()
        self.nb_bonds = len(self._bonds.id.values)
        self.timesteps = self._bonds.ts.values
        self.ids = self._bonds.id.values