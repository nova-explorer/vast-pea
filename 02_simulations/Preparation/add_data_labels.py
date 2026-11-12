#!/bin/python3

"""
Small custom utility to copy type labels from one file to another

Run using:
    ./add_data_labels.py <file without labels> <file with labels>
    ./add_data_labels.py system_equil.data system_cell.data

Outputs <file without labels>2
        system_equil.data2
"""

import sys
import os

def main(input_file, reference_file, path):
    header, body = read_input(input_file, path)

    labels = read_reference(reference_file, path)

    write_data(input_file, header, labels, body)

def read_input(input_file, path):
    with open(path+input_file, 'rt') as file_in:
        lines = file_in.readlines()

    for i, line in enumerate(lines):
        if "Masses" in line: # header end trigger
            header_end = i-1
            break
    return lines[:header_end], lines[header_end+1:]

def read_reference(reference_file, path):
    with open(path+reference_file, 'rt') as file_in:
        lines = file_in.readlines()

    for i, line in enumerate(lines):
        if "Atom Type Labels" in line: # start trigger
            start = i
        elif "Masses" in line:
            end = i-1
            break

    return lines[start:end]

def write_data(filename, header, labels, body):
    with open(path + filename + "2", "wt") as file_out:
        file_out.write("".join(header) + "\n")
        file_out.write("".join(labels) + "\n")
        file_out.write("".join(body) + "\n")

if __name__ == "__main__":
    input_file = sys.argv[1]
    template_file = sys.argv[2]
    path = os.getcwd() + "/"

    main(input_file, template_file, path)