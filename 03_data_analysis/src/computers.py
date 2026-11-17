# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: Functions to compute certain properties

TODO:
    -
"""

import numpy as np
from scipy.spatial.distance import cdist
import xarray as xr
from copy import copy

from constants import *

def compute_distance_matrix(positions1: np.ndarray, positions2: np.ndarray, bounds: np.ndarray) -> np.ndarray:

    nb_atoms = len(positions1)
    comps = len(positions1[0])

    distance_matrix = np.zeros([len(positions1), len(positions2)])

    for comp_id in range(comps):
        comp_distance = cdist( positions1[:, [comp_id]], positions2[:, [comp_id]], "euclidean" )

        distance_matrix += ( comp_distance - bounds[comp_id] * np.round(comp_distance/bounds[comp_id]) )**2

    if np.any(distance_matrix<0):
        raise ValueError("Negative values in distance matrix")

    if np.isnan(distance_matrix).any():
        raise ValueError("NaN values in distance matrix")

    if np.any(distance_matrix<1e-2):
        distance_matrix = np.where(distance_matrix<1e-2, 0, distance_matrix)

    return np.sqrt(distance_matrix)

def compute_rdf(pair_positions1: xr.DataArray, pair_positions2: xr.DataArray, bounds: xr.DataArray, cutoff: float = 25, n_bins: int = 20, low_cutoff: float = 0, print_flag=True) -> xr.DataArray:

    # check if pair_positions1 and pair_positions2 are compatible
    if len(pair_positions2.ts.values) != len(pair_positions1.ts.values):
        raise IndexError("Pairs do not have the same timesteps")
    if len(pair_positions1.comp.values) != len(pair_positions2.comp.values):
        raise IndexError("Pairs do not have the same components")
    if pair_positions1.isnull().any() or pair_positions2.isnull().any():
        raise Warning("Pair positions contain NaNs")

    increment = cutoff / n_bins
    bins = np.linspace(increment/2, cutoff-increment/2, n_bins)
    edges = np.linspace(0, cutoff, n_bins+1)

    hollow_sphere_volume = 4/3 * np.pi * ( (bins+increment/2)**3 - (bins-increment/2)**3 )


    total_ts = len(pair_positions1.ts.values)
    total_id_1 = len(pair_positions1.id.values)
    total_id_2 = len(pair_positions2.id.values)
    total_id = (total_id_1+total_id_2)/2

    rdf_traj = []
    for cnt_ts in range(total_ts):
        if print_flag:
            print('\r\tComputing RDF on timestep {}/{}'.format
                  (cnt_ts+1, total_ts),
                   end="\r")
        cell_volume = bounds.isel(ts=cnt_ts).prod().values
        normalize_factor = cell_volume / total_id**2 / hollow_sphere_volume

        rdf_array = np.zeros(n_bins)
        dist = compute_distance_matrix(pair_positions1.isel(ts=cnt_ts).transpose('id', 'comp').values,
                                       pair_positions2.isel(ts=cnt_ts).transpose('id', 'comp').values,
                                       bounds.isel(ts=cnt_ts).values
                                       )
        for cnt_id in range(total_id_1):

            data = dist[cnt_id]
            data = data[(data>low_cutoff) & (data<=cutoff)]

            bin_count, _ = np.histogram(data, edges)

            rdf_array += bin_count

        rdf_traj.append(rdf_array * normalize_factor)

    if print_flag:
        print("\n", end="\r")

    rdf = xr.DataArray(rdf_traj, coords = [pair_positions1.ts.values, list(bins)], dims = ['ts', 'distance'])

    return rdf

def compute_msd(positions: xr.DataArray):

    comps = len(positions.comp.values)
    nb_atoms = len(positions.id.values)
    ts = len(positions.ts.values)

    # bounds = bounds.transpose('ts', 'comp').values
    positions = positions.transpose('ts', 'id', 'comp').values # pyright: ignore[reportAssignmentType]

    displacement = positions - positions[0]

    total_displacement = np.zeros([ts, nb_atoms])

    for i in range(comps):
        total_displacement += ( displacement[:,:,i] )**2

    return total_displacement.mean(axis=1)

def compute_msd_cm(positions: xr.DataArray, types: list, masses: list) -> np.ndarray:
    """
    positions: ts, id, comp
        [                # ts
            [            # id
                [ x y z] # comp
            ]
        ]

    types: id
        [             # id
            type (int)
        ]

    masses: [m0, m1, m2, m3, ...] floats
    """

    comps = len(positions.comp.values)
    nb_atoms = len(positions.id.values)
    ts = len(positions.ts.values)
    # ts_list = positions.ts.values

    positions = positions.transpose('ts','id','comp').values # pyright: ignore[reportAssignmentType]

    mass_array = []
    for i in range(nb_atoms):
        mass_array.append(masses[types[i]-1])
    mass_array = np.array(mass_array)
    total_mass = mass_array.sum()

    cm = (positions * mass_array[None,:,None]).sum(axis=1) / total_mass

    cm_displacement = cm - cm[0]

    displacement = positions - positions[0] - cm_displacement[:,None,:]

    total_displacement = np.zeros([ts, nb_atoms])

    for i in range(comps):
        total_displacement += ( displacement[:,:,i] )**2

    return total_displacement.mean(axis=1)

def compute_bond_autocorrelation(bonds):

    b_ac = []
    bonds_intersect = bonds[0]
    nb_bonds = len(bonds[0])

    for ts in bonds:

        bonds_intersect = bonds_intersect.intersection(ts)
        b_ac.append(len(bonds_intersect)/nb_bonds)

    return b_ac

def compute_nb_reactions(bonds):
    nb_reactions = []
    nb_bonds = len(bonds[0])
    old_bonds = bonds[0]

    for ts in bonds:

        nb_reactions.append( nb_bonds - len(old_bonds.intersection(ts)) )
        old_bonds = ts

    return nb_reactions

def compute_crosslinking_density(lammps_data, effective=True, volume=None, bounds=None, crosslink_point_type=10, crosslink_types=[3,5], is_poss=False, poss_type=15, type_per_poss=8, poss_offset=108):
    if (not volume and not bounds) or (volume and bounds) :
        raise ValueError("Need to provide bounds OR volume")
    elif volume:
        vol = volume
    elif bounds:
        vol = np.prod([ bounds[i][1] - bounds[i][0] for i in bounds ])

    sub_atoms = lammps_data.get_atoms().where(
        (lammps_data.get_atoms().type == crosslink_point_type) |
        (lammps_data.get_atoms().type.isin(crosslink_types)),
        drop=True
    )
    sub_data = copy(lammps_data)
    sub_data.set_atoms(sub_atoms)
    sub_data.establish_connections_subset()

    crosslink_point_ids = np.sort(sub_data.get_atoms().where(sub_data.get_atoms().type==crosslink_point_type, drop=True).id.values)

    if is_poss:
        poss_ids = np.sort(lammps_data.get_atoms().where(lammps_data.get_atoms().type==poss_type, drop=True).id.values) [::type_per_poss]
        poss_range_ids = []
        for i in poss_ids:
            poss_range_ids += list(range(i, i+poss_offset))

        hardener_point_ids = np.array([i for i in crosslink_point_ids if i not in poss_range_ids])
        hardener_point_ids = hardener_point_ids.reshape(int(len(hardener_point_ids)/2),2)

        poss_point_ids = np.array([i for i in crosslink_point_ids if i in poss_range_ids])

    else:
        hardener_point_ids = crosslink_point_ids.reshape(int(len(crosslink_point_ids)/2),2)

    def count_crosslinks_hardener(a):
        return len(sub_data.bonds[a[0]]) + len(sub_data.bonds[a[1]])
    def count_crosslinks_poss(a):
        return len(sub_data.bonds[a])

    def count_effective_hardener(a):
        # 1 for AFD + 0.5 for each amine crosslink
        if a < 3:
            return 0 # not elastically effective if nb of crosslink on AFD < 3
        elif a == 3:
            return 2.5
        elif a == 4:
            return 3

    def count_non_effective_hardener(a):
        # 1 for each amine with 2 crosslinks
        if a < 2:
            return 0
        elif a < 4:
            return 1
        elif a == 4:
            return 2

    def count_effective_poss(a):
        # 0.5*3 for POSS corner + 1 for benzene + 0.5 for each amine crosslink
        if a == 0:
            return 0 # not elastically effective if not grafted to network
        elif a == 1:
            return 3
        elif a == 2:
            return 3.5

    def count_non_effective_poss(a):
        # 1 for each amine with 2 crosslinks
        if a < 2:
            return 0
        elif a == 2:
            return 1

    count_hard = np.apply_along_axis(count_crosslinks_hardener, 1, hardener_point_ids)
    if effective:
        count_hard = list(map(count_effective_hardener, count_hard))
    else:
        count_hard = list(map(count_non_effective_hardener, count_hard))
    count = sum(count_hard)

    if is_poss:
        count_poss = list(map(count_crosslinks_poss, poss_point_ids))
        if effective:
            count_poss = list(map(count_effective_poss, count_poss))
        else:
            count_poss = list(map(count_non_effective_poss, count_poss))
        count += sum(count_poss)

    return count / vol / ANGSTROM_TO_METER**3 / AVOGADRO_NUMBER

def compute_legendre_poly(vectors:xr.DataArray, order:int, origin_nb:int) -> np.ndarray:

    if len(vectors.ts.values) < origin_nb:
        raise IndexError("Nb of origins is larger than number of data")

    if origin_nb > 1:
        origins = sorted(sample(list(vectors.ts.values), origin_nb))
    elif origin_nb == 1:
        origins = [0]
    else:
        raise IndexError("Nb of origins not valid")

    oacf = []
    for origin in origins:

        vectors_o = vectors.where(vectors.ts >= origin)

        dot_prod = vectors_o.dot( vectors_o.isel(ts=0), dim='comp' )
        if order == 1:
            oacf.append( dot_prod.mean(dim='id') )

        elif order == 2:
            oacf.append( (dot_prod**2).mean(dim='id') )

        else:
            raise IndexError("Legendre polynomial order not implemented")

    max_len = len(oacf[0])
    oacf = np.nanmean( [ np.append(i, [np.nan]*(max_len-len(i))) for i in oacf ], axis=0 )

    if order == 2:
        oacf = 3/2 * oacf - 1/2

    return oacf


def compute_skewness_voxel(coms, bounds, vox_per_dim: int, verbose=True):
    timesteps = coms.ts.values
    ids = coms.id.values
    comps = coms.comp.values

    bounds = bounds.mean(dim='comp')
    coms = coms.positions
    coms -= bounds * np.floor(coms/bounds)

    n_voxels = vox_per_dim ** 3
    vox_half_sizes = bounds / 2 / vox_per_dim

    if verbose:
        cell_len = bounds
        cell_vol = bounds**3
        vox_len = 2*vox_half_sizes
        vox_vol = (2*vox_half_sizes)**3
        print("Simulation cell volume: {:.2f}±{:.2f}, length: {:.2f}±{:.2f}".format(
            cell_vol.mean(dim='ts'),
            cell_vol.std( dim='ts'),
            cell_len.mean(dim='ts'),
            cell_len.std( dim='ts'),
        ))
        print("Voxel volume: {:.2f}±{:.2f}, length: {:.2f}±{:.2f}".format(
            vox_vol.mean(dim='ts'),
            vox_vol.std( dim='ts'),
            vox_len.mean(dim='ts'),
            vox_len.std( dim='ts'),
        ))
        print("Simulation cell - voxels volume: {:.5f}±{:.5f}, length: {:.5f}±{:.5f}".format(
            (cell_vol - n_voxels*vox_vol).mean(dim='ts'),
            (cell_vol - n_voxels*vox_vol).std( dim='ts'),
            (cell_len - vox_per_dim*vox_len).mean(dim='ts'),
            (cell_len - vox_per_dim*vox_len).std( dim='ts'),
        ))

    voxels = xr.Dataset()
    voxels['centers'] = xr.DataArray(
        np.array([ [ [ ( (2*i+1)*vox_half_sizes, (2*j+1)*vox_half_sizes, (2*k+1)*vox_half_sizes) for k in range(vox_per_dim)] for j in range(vox_per_dim)] for i in range(vox_per_dim)]).reshape(n_voxels, 3, len(timesteps)),
        coords=[range(n_voxels), comps, timesteps],
        dims=["id", "comp", "ts"]
        )

    count_ts = []

    vox_half_sizes_arr = vox_half_sizes.transpose('ts').values
    coms_arr = coms.transpose('ts','id','comp').values
    centers_arr = voxels.centers.transpose('ts', 'id', 'comp').values
    for ts in range(len(timesteps)):

        count = [0] * len(voxels.id)
        for vox_id in range(len(voxels.id)):

            for poss_id in range(len(ids)):

                count[int(vox_id)] += is_point_in_voxel(
                    point_center=coms_arr[ts, poss_id, :],
                    voxel_center=centers_arr[ts, vox_id, :],
                    voxel_half_length=vox_half_sizes_arr[ts]
                    )

        count_ts.append(count)
    voxels['nq'] = xr.DataArray(count_ts, coords=[ voxels.ts,voxels.id], dims=["ts", "id"])

    mean_nq = voxels.nq.mean(dim='id')
    std_nq = voxels.nq.std(dim='id')
    skewness = n_voxels / (n_voxels-1) / (n_voxels-2) * ( ((voxels.nq - mean_nq)/std_nq)**3 ).sum(dim='id')

    max_q = int(voxels['nq'].max().values)
    q_count = [ voxels['nq'].where(voxels['nq']==i).count(dim='id') for i in range(max_q+1) ]

    return skewness, xr.concat(q_count, dim='q')

def is_point_in_voxel(point_center:np.ndarray, voxel_center: np.ndarray, voxel_half_length: float) -> bool:
    dx = abs(point_center[0] - voxel_center[0]) - voxel_half_length
    dy = abs(point_center[1] - voxel_center[1]) - voxel_half_length
    dz = abs(point_center[2] - voxel_center[2]) - voxel_half_length

    return np.max([dx, dy, dz]) < 0