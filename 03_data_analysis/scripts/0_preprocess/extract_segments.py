from vast_pea.lammps_traj import trajectory
from vast_pea.lammps_segments import segment
from vast_pea.masses import PCFF_MASSES, types_to_masses

traj = trajectory(path="../../data/BER/save/",
                  pattern="afd_325K_all.traj_data.atoms.nc", # we use a saved netCDF file since faster and smaller
                  filter_type=False, # when restoring, filter doesn't have an effect
                  filter=False,
                  restore=True,
                  verbose=True
                  )

# a list of the types in the system in the order of their atom type id
types = ['c',   # 1
         'c1',  # 2
         'c2',  # 3
         'c3',  # 4
         'cp',  # 5
         'h',   # 6
         'hc',  # 7
         'hn',  # 8
         'ho2', # 9
         'nn',  # 10
         'oc',  # 11
         'oh',  # 12
         's1',  # 13
        ]
# transfer that list of types to a list of masses
masses = types_to_masses(types, PCFF_MASSES)

# define the segments
segment_identifiers = [
    { # AFD 1st segment
        # types of the identifiable beginning and end of the segment
        'end_types':[10,13], # here identified beginnings and ends of AFD segments are nn and s1
        # offset to the end_types for the actual beginning and end of the segment
        'end_offset':[0,0], # no offset to the atom ids of identified beginnings and ends is needed for this segment
        # atom ids to include in the center of mass calculation of the segment relative to the beginning of the segment
        # A good way to find this is checking the atom ids we want for each segment is checking the atom ids of monomers in data or mol2 files
        'com_range':["-5:6"], # this here means that CoM for AFD segments consider from atom ids beginning-5 up to beginning+6
        # offset in atom id to not scan for a new segment after the end of the previous segment
        'exclusivity_offset':False # not needed here
    },
    { # AFD 2nd segment
        'end_types': [13,10], 'end_offset':[0,0], 'com_range':["0:11"], 'exclusivity_offset':False
    },
    { # DGEBA 1rst segment
        'end_types':[11,1], # here identified beginning and ends of DGEBA segments are c and oc
        'end_offset':[0,0], # no offset is needed
        'com_range':["0:7","24:27"], # CoM for DGEBA segments consider from atom ids beginning up to beginning+7, and from beginning+24 up to beginning+27
        'exclusivity_offset':False # not needed here
    },
    { # DGEBA/2nd segment
        'end_types': [1,11], 'end_offset':[0,0], 'com_range':["0:7","23:26"], 'exclusivity_offset':False
    }
    ]

# calculate segments
segment_traj = segment(traj=traj, # trajectory object
                       mass_list=masses, # list of masses for every type
                       identifiers=segment_identifiers, # list of dicts detailing segments
                       verbose=True # print progress on calculating segemnets
                       )

# output the center of mass of segments and all atoms in the trajectory to a LAMMPS dump file.
# It is fairly slow, so it only outputs the first timestep.
# segment_traj.dump_unwrapped_coms_and_atoms(path="../../data/BER/",
#                                          name="afd_ber_325K_segments.dump",
#                                          traj=traj
#                                          )

# this will save 4 files
#   - <name>.segment_data.coms.nc: xr.Dataset with center of masses of segments, saved in netCDF format
#   - <name>.segment_data.details.dc: printed identifiers argument
#   - <name>.segment_data.ends.nc: xr.Dataset with information about the beginnings and ends of segments, saved in netCDF format
#   - <name>.segment_data.vectors.nc: xr.Dataset with beginning->end vectors of segments, saved in netCDF format
segment_traj.save(name="afd_325K",
                  path="../../data/BER/save/"
                  )