#!/usr/bin/env python

from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

parser = argparse.ArgumentParser(description="given a list of silva fasta files pull out the taxids")
parser.add_argument('fasta', help='path to fasta file containing silva sequences to be matched')
args = parser.parse_args()

# Open the fasta database
silva_fa = open(args.fasta, "r")
# Open the untrimmed taxa mapping file
Silva_tax_map = open("Silva_108_taxa_mapping.txt", "r")
# Write a file to store the trimmed taxa mapping data
outname = args.fasta.rstrip(".fasta") + ".txt"
Silva_tax_map_trimmed = open(outname, "w")

# Create an empty dictionary of sequence names
name_dict = {}

# For each line in the taxa mapping file
for line in Silva_tax_map:
    # get the name and store it as a variable
    namestop = line.find("\t")
    name = line[:namestop]
    # get the taxonomic info and store it as a variable
    tax = line[namestop:]
    # store the name and the taxonomic assignment as a dictionary item
    f = name, tax
    name_dict[(name, tax)] = f

# For each sequence in the trimmed fasta file
for title, seq in SimpleFastaParser(silva_fa):
    # Iterate through the dictionary of ID:taxonomic_info
    for (name, tax) in name_dict:
        # if the sequnce title matches the name in the dictionary...
        if name == title:
            # write the sequence ID and tax record to the trimmed taxa file
            name = name.rstrip("\n")
            Silva_tax_map_trimmed.write("%s%s\n" % (name, tax))
