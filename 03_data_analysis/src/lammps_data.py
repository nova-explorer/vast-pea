import xarray as xr
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist


import io_

class data_file():
    """
    TODO:
        - __read_data keep only numbers could be handled more dynamically. LAMMPS usually has a pretty rigid formatting for sections spacing tho.
    """
    def __init__(self,
                 path="./", pattern="*.data", verbose=True
                 ) -> None:
        self.verbose = verbose

        file = io_.create_file_list(path, pattern)
        if len(file) > 1:
            # TODO clean
            raise ValueError("More than 1 file matching pattern.")
        self.file = str(file[0])

        self.__print("\tReading data file....\n")
        self._velocities = False
        self.__read_data()

        """
        self._header list of str
        self._labels list of str
        self._masses tuple
        self._coeffs list of str
        self._atoms xarray for compatibility
        self._velocities list of str
        self._interactions list of str
        """

        self.nb_atoms = len(self._atoms.id.values)
        self.ids = self._atoms.id.values
        self.comps = ['x', 'y', 'z']

        self.__add_masses()

        self.__print("\tData file done!\n")

    def __read_data(self) -> None:

        # Reading the file

        lines = io_.read_file(self.file)

        header_start = 0
        bounds_start = 0
        labels_start = 0
        masses_start = 0
        coeffs_start = 0
        atoms_start = 0
        interactions_start = 0

        velocities_start = False

        for i, line in enumerate(lines):
            if "xlo xhi" in line:
                bounds_start = i
            elif "Atom Type Labels" in line:
                labels_start = i
            elif "Masses" in line:
                masses_start = i
            elif "Pair Coeffs" in line:
                coeffs_start = i
            elif "Atoms # full" in line:
                atoms_start = i
            elif "Velocities" in line:
                velocities_start = i
            elif "Bonds" in line:
                interactions_start = i

        self._header = lines[header_start:bounds_start]
        self._labels = lines[labels_start:masses_start]
        self._coeffs = lines[coeffs_start:atoms_start]
        self._interactions = lines[interactions_start:]

        if velocities_start:
            self._velocities = lines [velocities_start:interactions_start]

        # Bounds to dict(list)

        bounds = lines[bounds_start:labels_start-1]
        self._bounds = { i.rstrip().split()[2][0]: [ float(i.rstrip().split()[0]), float(i.rstrip().split()[1]) ] for i in bounds}

        # Masses to dict
        masses = lines[masses_start+2:coeffs_start-1] # only keep numbers
        if '#' in masses[0]:
            self._masses = {i.rstrip().split()[-1]: float(i.rstrip().split()[1]) for i in masses }
        else:
            self._masses = {float(i.rstrip().split()[0]): float(i.rstrip().split()[1]) for i in masses }

        # Atoms to xarray
        if velocities_start:
            atoms = lines[atoms_start+2:velocities_start-1] # only keep numbers
        else:
            atoms = lines[atoms_start+2:interactions_start-1] # only keep numbers

        df = self.__lines_to_df(atoms)

        self._atoms = df.to_xarray()

    def __lines_to_df(self, lines) -> pd.DataFrame:

        columns = ['id', 'molecule', 'type', 'charge', 'x', 'y', 'z', 'ix', 'iy', 'iz'] # from style full

        lines = [ i.rstrip().split()[:-2] if '#' in i else i.rstrip().split() for i in lines ] # remove comments at position -2

        df = pd.DataFrame(lines, columns=columns).apply(pd.to_numeric)
        df = df.set_index('id')

        return df

    def __add_masses(self) -> None:
        mass_list = list(self._masses.values())
        masses = []

        for i in range(self.nb_atoms):
            type = int(self._atoms['type'].isel(id=i))
            masses.append(mass_list[type-1])

        self._atoms['mass'] = xr.DataArray(masses, coords=[self.ids], dims=['id'])

    def __extract_bonds(self):
        start = 0
        end = 0

        for i, line in enumerate(self._interactions):
            if "Bonds" in line:
                start = i + 2
            elif "Angles" in line:
                end = i - 1
                break

        return self._interactions[start:end]

    def __print(self, text) -> None:
        if self.verbose:
            print(text, end="\r")

    def unwrap(self) -> None:
        for i in self.comps:
            self._atoms[i] += self._atoms['i'+i] * (self._bounds[i][1] - self._bounds[i][0])
            self._atoms['i'+i] *= 0

    def agglomerate(self) -> None:

        positions = self.get_positions().to_numpy()

        for comp_id in range(len(self.comps)):
            comp = self.comps[comp_id]

            # dist = cdist(positions[:,[comp_id]], positions[:,[comp_id]], "euclidean") [0]
            dist = positions[:,comp_id] - [positions[0,comp_id]]

            bound = self._bounds[comp][1] - self._bounds[comp][0]

            self._atoms[comp] -= bound * np.round(dist/bound)

    def agglomerate_molid(self) -> None:

        for comp_id, comp in enumerate(self.comps):

            new_pos = xr.DataArray([],coords=[[]],dims=['id'])

            for molid in range(int(self._atoms.molecule.max())):

                if self._atoms.where(self._atoms.molecule == molid+1, drop=True).id.size > 0:

                    positions = self._atoms.where(self._atoms.molecule==molid+1, drop=True)[comp].to_numpy()

                    dist = positions - positions[0]

                    bound = self._bounds[comp][1] - self._bounds[comp][0]
                    positions = xr.DataArray(positions - bound * np.round(dist/bound), coords=[self._atoms.where(self._atoms.molecule==molid+1, drop=True).id], dims=['id'])

                    new_pos = xr.concat([new_pos, positions], dim='id')

            self._atoms[comp] = new_pos

    def recenter(self) -> None:
        for i in self.comps:
            com = (self._atoms[i] * self._atoms['mass']).sum() / self._atoms['mass'].sum()
            self._atoms[i] -= com

    def reset_bounds(self) -> None:
        for i in self.comps:
            bound = self._bounds[i][1] - self._bounds[i][0]
            self._bounds[i] = [-bound/2, bound/2]

    def reset_images(self) -> None:
        for i in self.comps:
            self._atoms['i'+i] *= 0

    def shrink(self) -> None:
        outmost = float(abs(self.get_positions()).max()) + 1 # pyright: ignore[reportArgumentType]

        for comp in self.comps:
            self._bounds[comp] = [-outmost, outmost]

    def write_file(self, filename) -> None:
        io_.to_data(filename, self._header, self._bounds, self._labels, self._masses, self._coeffs, self._atoms, self._velocities, self._interactions)

    def establish_connections(self):
        bond_lines = self.__extract_bonds()

        bond_data = {i:[] for i in self._atoms.id.values}
        for line in bond_lines:

            atom1, atom2 = [int(i) for i in line.rstrip().split(" ") [2:4]]
            bond_data[atom1].append(atom2)
            bond_data[atom2].append(atom1)

        self.bonds = bond_data

    def establish_connections_subset(self):
        bond_lines = self.__extract_bonds()

        atom_ids = self._atoms.id.values

        bond_data = {i:[] for i in atom_ids}
        for line in bond_lines:

            atom1, atom2 = [int(i) for i in line.rstrip().split(" ") [2:4]]

            if atom1 in atom_ids and atom2 in atom_ids:
                bond_data[atom1].append(atom2)
                bond_data[atom2].append(atom1)

        self.bonds = bond_data

    def clusterize(self, atoms):

        self.avail_atoms = list(atoms.id.values)
        self.clusters = []

        for i in atoms.id.values:
            current_cluster = []
            current_cluster = self.__connect(i, current_cluster)

            if current_cluster:
                self.clusters.append(current_cluster)

        current_mol = int(self._atoms.molecule.max()) + 1
        clustered_frame = atoms.to_dataframe()

        for clust_i in self.clusters:
            for atom_j in clust_i:
                clustered_frame.loc[atom_j, 'molecule'] = current_mol
            current_mol += 1

        return clustered_frame.to_xarray()

    def __connect(self, atom_id, current_cluster):
        if atom_id in self.avail_atoms:

            self.avail_atoms.remove(atom_id)
            current_cluster.append(atom_id)

            for connection in self.bonds[atom_id]:
                current_cluster = self.__connect(connection, current_cluster)

        return current_cluster

    def get_header(self) -> list:
        return self._header
    def get_bounds(self) -> dict:
        return self._bounds
    def get_labels(self) -> list:
        return self._labels
    def get_masses(self) -> dict:
        return self._masses
    def get_coeffs(self) -> list:
        return self._coeffs
    def get_atoms(self) -> xr.Dataset:
        return self._atoms
    def get_interactions(self) -> list:
        return self._interactions

    def get_positions(self) -> xr.DataArray:
        return xr.DataArray(np.transpose([self._atoms.x, self._atoms.y, self._atoms.z]), coords=[self.ids, self.comps], dims=['id', 'comp'])

    def set_atoms(self, new_atoms) -> None:
        self._atoms = new_atoms.copy()
        self._atoms['molecule'] = self._atoms['molecule'].astype(int)
        self._atoms['type'] = self._atoms['type'].astype(int)
        self._atoms['ix'] = self._atoms['ix'].astype(int)
        self._atoms['iy'] = self._atoms['iy'].astype(int)
        self._atoms['iz'] = self._atoms['iz'].astype(int)

        self.nb_atoms = len(self._atoms.id.values)
        self.ids = self._atoms.id.values
