#!/bin/python3

"""
Small custom utility to rescale the charges of specified atoms so the net charge of the system is neutral.
The hc (hydrogen bonded to carbon) type is typically chosen since there's a lot of those and the charge offset will be minimal.
More than one type could be specified.

Run using:
    ./charge_rescale.py <file> <atom types to rescale>
    ./charge_rescale.py system_poly_90.data 7

Outputs <file rescaled>2
        system_poly_90.data2
"""

import sys
import os

def read_file(filename):
    with open(filename, 'rt') as file_in:
        lines = file_in.readlines()

    for i, line in enumerate(lines):
        if "Atoms # full" in line:
            start = i+2
        elif "Velocities" in line or "Bonds" in line:
            end = i-1
            break
    head = lines[:start]
    data = lines[start:end]
    tail = lines[end:]
    return head, data, tail

def format_data(data):
    data = [i.split(" ") for i in data]
    data = [ [ convert(j, val) for j, val in enumerate(i) ] for i in data ]
    return data

def convert(pos, val):
    if pos in [0, 1, 2, 7, 8, 9]:
        return int(val)
    elif pos in [3, 4, 5, 6]:
        return float(val)

def get_net_charge(data):
    data = [ i[3] for i in data ]
    charge = sum(data)
    return charge

def get_nb_typec(data, typec):
    data = [ i for i in data if i[2] in typec ]
    nb_typec = len(data)
    return nb_typec

def apply_rescale(data, typec, rescale):
    data = [ [ j+rescale if j==i[3] else j for j in i ] if i[2] in typec else i for i in data ]
    return data

def print_to_file(head, data, tail, filename, path):
    filename_out = path + filename + "2"

    data = [ [str(j) for j in i] for i in data ]
    data = [ " ".join(i) + "\n" for i in data ]
    print(filename_out)
    with open(filename_out, 'wt') as file_out:
        file_out.write("".join(head))
        file_out.write("".join(data))
        file_out.write("".join(tail))

def check_new_file(file, path, typec):
    file = file.split(".")[0] + "_2." + file.split(".")[1]
    head, data, tail = read_file(path+file)
    data = format_data(data)

    charge = get_net_charge(data)
    nb_typec = get_nb_typec(data, typec)
    rescale = charge / nb_typec
    print("New file")
    print("\tNet charge: {}".format(charge))
    print("\tNb of type {} atoms: {}".format(typec, nb_typec))
    print("\tRescaling offset: {}".format(rescale))

def main(file, path, typec):
    head, data, tail = read_file(path+file)
    data = format_data(data)

    charge = get_net_charge(data)
    nb_typec = get_nb_typec(data, typec)
    rescale = -charge / nb_typec
    print("Current file")
    print("\tNet charge: {}".format(charge))
    print("\tNb of type {} atoms: {}".format(typec, nb_typec))
    print("\tRescaling offset: {}".format(rescale))

    data = apply_rescale(data, typec, rescale)

    print_to_file(head, data, tail, file, path)

    check_new_file(file, path, typec)

if __name__ == "__main__":
    input_file = sys.argv[1]
    path = os.getcwd()+"/"
    typec = sys.argv[2:]
    typec = [int(i) for i in typec]

    main(input_file, path, typec)
