#!/usr/bin/env python
# Convert fastq files to fasta.

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='convert fastq to fasta')
parser.add_argument('fastq', help='path to fasta file to be converted')
args = parser.parse_args()

fastq = str(args.fastq)
name = fastq.replace(".fastq", "")

fasta = SeqIO.convert(args.fastq, "fastq", name + ".fasta", "fasta")
