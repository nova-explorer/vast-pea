import xarray as xr
import numpy as np
import itertools
import os

import io_

class segment():
    def __init__(self,
                traj,
                mass_list:list,
                identifiers:list,
                restore=False,
                restore_path=None,
                restore_pattern=None,
                verbose=True
                 ) -> None:

        if verbose: print("\tCalculating segments!")
        if not restore:
            self.identifiers = identifiers
            self.comps = traj.comps
            self.timesteps = traj.timesteps

            self._ends = self.__make_id_map(traj.get_atoms().isel(ts=0))

            self.ids = self._ends.id.values
            self.nb_segments = len(self.ids)
            self.cases = range(int(self._ends.case.max()) + 1)

            self._vectors = self.__compute_vectors(traj)
            self._coms = self.__compute_centers_of_mass(traj, mass_list)

        else:
            file_list = io_.create_file_list(restore_path, restore_pattern)
            self.__restore(file_list)

            self.comps = self._coms.comp
            self.timesteps = self._coms.ts.values
            self.ids = self._ends.id.values
            self.nb_segments = len(self.ids)
            self.cases = range(int(self._ends.case.max()) + 1)

    def __make_id_map(self, atom_frame) -> xr.Dataset:

        pairs = self.__find_pairs(atom_frame)
        flattened_pairs = [ pair for case in pairs for pair in case ]
        flattened_cases = [ case_nb for case_nb, case in enumerate(pairs) for _ in case ]

        end_ids = xr.Dataset(
            {"pairs": (["id", "pair_nb"], flattened_pairs),
             "case": (["id"], flattened_cases)
            },
            coords={"id": range(len(flattened_cases)),
                    "pair_nb": range(2),
                    }
        )

        end_ids['end'] = end_ids.pairs + [ self.identifiers[i]['end_offset'] for i in end_ids.case.values ]

        return end_ids

    def __find_pairs(self, atoms) -> np.ndarray:
        atom_types = atoms.type

        all_pairs = []
        for identifier in self.identifiers:
            matching_types = list(set(identifier['end_types']))
            matching_atoms = atom_types.where(atom_types.isin(matching_types), drop=True).id.values

            pairs = [ [matching_atoms[i], matching_atoms[i+1]] for i in range(len(matching_atoms)-1) ]
            pairs = [ pair for pair in pairs if list(atom_types.sel(id=pair).values) == identifier['end_types'] ]

            if identifier['exclusivity_offset']:
                chosen = [pairs[0]]
                for i in pairs:
                    if i[0] > chosen[-1][0]+identifier['exclusivity_offset']:
                        chosen.append(i)
                pairs = chosen
            all_pairs.append(pairs)

        return all_pairs

    def __compute_vectors(self, traj) -> xr.Dataset:
        positions = traj.get_positions()
        bounds = traj.get_bounds()

        atom0 = positions.sel(id=self._ends.end.sel(pair_nb=0)).drop_vars(names=['pair_nb'])
        atom1 = positions.sel(id=self._ends.end.sel(pair_nb=1)).drop_vars(names=['pair_nb'])

        atom0['id'] = range(len(atom0))
        atom1['id'] = range(len(atom1))

        vector = atom1 - atom0

        vector = vector - bounds * np.round(vector/bounds)

        length = np.sqrt(vector.dot(vector, dim='comp'))
        vector /= length

        ds = xr.Dataset(
            data_vars = {
                'vector': vector,
                'length': length,
                'pos0': atom0,
                'pos1': atom1
            },
            coords = {
                'id': atom0.id.values,
                'ts': traj.timesteps,
                'comp': traj.comps
            }
        )

        return ds

    def __compute_centers_of_mass(self, traj, mass_list) -> xr.Dataset:

        com_traj = []
        for i, identifier in enumerate(self.identifiers):
            segment_traj = self.__get_segment_traj(traj.get_unwrapped(), identifier['com_range'], i)

            segment_traj['group_mass'] = self.__get_mass_group(
                traj.get_atoms().type.isel(ts=0).values.astype(int),
                mass_list,
                segment_traj.group_ids,
                list(traj.get_atoms().id.values)
                )

            segment_traj['mass'] = segment_traj['group_mass'].sum(dim='id_group')

            segment_traj['positions'] = (segment_traj['group_positions'] * segment_traj['group_mass']).sum(dim='id_group') / segment_traj['mass']

            segment_traj['type'] = xr.DataArray(
                [i] * len(segment_traj.id.values),
                coords = [segment_traj.id.values],
                dims = ['id']
            )

            segment_traj = segment_traj.drop_vars(["group_positions", "group_ids", "group_mass"]).drop_dims(['id_group'])
            com_traj.append(segment_traj)

        com_traj = xr.concat(com_traj, dim="id", coords="different", compat='equals' )

        com_traj['displacement'] = com_traj['positions'] - com_traj['positions'].isel(ts=0)
        com_traj['displacement'] -= com_traj['positions'].mean(dim='id') - com_traj['positions'].isel(ts=0).mean(dim='id')

        return com_traj

    def __get_segment_traj(self, positions, com_steps, current_case) -> xr.Dataset:
        com_ranges = self.__get_com_ranges(com_steps, current_case)

        segment_traj = []
        segment_ids = []
        for range_i in com_ranges:
            segment_traj.append(positions.sel(id=range_i).values)
            segment_ids.append(list(range_i))

        segment_traj = xr.Dataset(
            {
                'group_positions': (["id", "id_group", "ts", "comp"], segment_traj),
                'group_ids'      : (["id", "id_group"], segment_ids)
            },
            coords = {
                "id": self._ends.where(self._ends.case==current_case, drop=True).id,
                "id_group": range(np.array(segment_traj).shape[1]),
                "ts": self.timesteps,
                "comp": self.comps
            }
        )
        return segment_traj

    def __get_com_ranges(self, com_steps, current_case) -> list:

        com_ranges = []
        for step in com_steps:

            if step.removeprefix('-').isnumeric():
                com_ranges.append(int(step))

            elif ':' in step and len(step.split(':')) == 2:
                step_range = [ int(i) for i in step.split(':') ]
                com_ranges.extend(range(step_range[0], step_range[1]+1))

            else:
                raise ValueError("com_range step not valid:", step)

        com_ids = self._ends.where(self._ends.case==current_case, drop=True).pairs.sel(pair_nb=0).values
        com_ranges = [ [i]*com_ids.shape[0] for i in com_ranges ]
        com_ids = (com_ids + com_ranges).T.astype(int)

        return com_ids

    def __get_mass_group(self, types, mass_list, group_ids, ids_index:list) -> xr.DataArray:
        mass_array = []
        for i in range(len(types)):
            mass_array.append(mass_list[types[i]-1])
        mass_array = np.array(mass_array)

        mass_group = []
        for i in group_ids.id.values:
            mass_group.append(
                [ mass_array[ids_index.index(j)] for j in group_ids.sel(id=i).values ]
            )

        mass_group = xr.DataArray(
            mass_group,
            coords = [group_ids.id.values, group_ids.id_group.values],
            dims = ['id', 'id_group']
        )
        return mass_group

    def __restore(self, file_list) -> None:
        CURRENT_MAX_FILES = 3
        currently_imported = {'ends':False, 'vectors':False, 'coms':False}

        if len(file_list) > CURRENT_MAX_FILES:
            raise EnvironmentError("Too many files for trajectory restore: "+", ".join(file_list))

        for filename in file_list:
            if ".segment_data.ends.nc" in filename:
                if not currently_imported['ends']:
                    self._ends = io_.read_dataset(filename)
                    currently_imported['ends'] = filename
                else:
                    self.__raise_imported("ends", currently_imported['ends'], filename)

            elif ".segment_data.vectors.nc" in filename:
                if not currently_imported['vectors']:
                    self._vectors = io_.read_dataset(filename)
                    currently_imported['vectors'] = filename
                else:
                    self.__raise_imported("vectors", currently_imported['vectors'], filename)

            elif ".segment_data.coms.nc" in filename:
                if not currently_imported['coms']:
                    self._coms = io_.read_dataset(filename)
                    currently_imported['coms'] = filename
                else:
                    self.__raise_imported("coms", currently_imported['coms'], filename)

            else:
                raise EnvironmentError("Filename doesn't match anything: ", filename)

    def __raise_imported(self, var_name, current_file, new_file) -> EnvironmentError:
        raise EnvironmentError(
            "Already imported {}: {}, current file: {}".format(var_name, current_file, new_file)
                )

    def save(self, name, path) -> None:
        filename_out = name
        details = self.identifiers

        io_.save_xarray(self._ends,
                       path,
                       filename_out+".segment_data.ends.nc"
                       )
        io_.save_xarray(self._vectors,
                       path,
                       filename_out+".segment_data.vectors.nc"
                       )
        io_.save_xarray(self._coms,
                       path,
                       filename_out+".segment_data.coms.nc"
                       )
        io_.save_dict(details,
                      path,
                      filename_out+".segment_data.details.dc")

    def dump_unwrapped_coms_and_atoms(self, path, name, traj) -> None:

        filename = os.path.join(path, name)
        atoms_template = "{id} 1 {type} {x} {y} {z} 0 0 0\n"
        cms_template = "{id} 1 {type} {x} {y} {z} 0 0 0\n"

        atoms = traj.get_atoms()
        positions = traj.get_unwrapped()

        cms_type_start = atoms.type.max().values+1
        n_all = len(atoms.id.values) + len(self._coms.id.values)
        n_atoms = len(atoms.id.values)
        with open(filename, 'wt') as f:
            for ts in [atoms.ts.values[0]]:
                f.write( 'ITEM: TIMESTEP\n' )
                f.write( str(ts)+'\n' )

                f.write( 'ITEM: NUMBER OF ATOMS\n' )
                f.write( str(n_all) + '\n' )

                f.write( 'ITEM: BOX BOUNDS pp pp pp\n' )
                for comp in atoms.bounds.comp.values:
                    f.write(
                        '0 ' + str(atoms.bounds.sel(ts=ts, comp=comp).values) + '\n'
                        )

                f.write( 'ITEM: ATOMS id mol type x y z ix iy iz\n')
                for i in atoms.id.values:
                    f.write(
                        atoms_template.format(
                            id=i,
                            type=atoms.sel(ts=ts, id=i).type.values,
                            x=positions.sel(ts=ts, id=i, comp='x').values,
                            y=positions.sel(ts=ts, id=i, comp='y').values,
                            z=positions.sel(ts=ts, id=i, comp='z').values,
                            )
                        )

                for i, id_i in enumerate(self._coms.id.values):
                    f.write(
                        cms_template.format(
                            id=i + n_atoms + 1,
                            type=cms_type_start + self._coms.type.sel(id=id_i).values,
                            x=self._coms.positions.sel(ts=ts, id=id_i, comp='x').values,
                            y=self._coms.positions.sel(ts=ts, id=id_i, comp='y').values,
                            z=self._coms.positions.sel(ts=ts, id=id_i, comp='z').values
                            )
                        )

    def get_coms(self):
        return self._coms
    def get_vectors(self):
        return self._vectors
    def get_ends(self):
        return self._ends