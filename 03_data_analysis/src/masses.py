
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

def types_to_masses(type_list, mass_dict) -> list:

    return [mass_dict[i] for i in type_list]