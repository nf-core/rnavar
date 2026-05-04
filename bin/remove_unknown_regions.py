#!/usr/bin/env python

import argparse
import os

def load_sequences_from_dict(dict_file):
    sequences = set()
    with open(dict_file, 'r') as file:
        for line in file:
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('SN:'):
                        sequences.add(part.split(':')[1])
    return sequences

def filter_bed_file(bed_file, sequences, output_file):
    with open(bed_file, 'r') as file, open(output_file, 'w') as out:
        for line in file:
            sequence = line.split('\t')[0]
            if sequence in sequences:
                out.write(line)

def main(bed_file, dict_file, output_file):
    sequences = load_sequences_from_dict(dict_file)
    filter_bed_file(bed_file, sequences, output_file)
    print(f"Output file {output_file} created in {os.getcwd()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove unknown regions from a BED file")
    parser.add_argument('bed', metavar='FILE', type=str, help="The BED file to remove unknown regions of")
    parser.add_argument('dict', metavar='FILE', type=str, help="The DICT file to get known regions from")
    parser.add_argument('output', metavar='FILE', type=str, help="The name of the output file")

    args = parser.parse_args()
    main(args.bed, args.dict, args.output)
