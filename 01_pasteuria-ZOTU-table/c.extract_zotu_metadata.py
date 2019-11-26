# Collect some sequence metadata which we'll use in data filtering

# Import modules.
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Set the file handle.
fasta_handle = open("pasteuria_alpha1_zotus_rerun.fasta")

# Open a blank output file. We're going to make a csv table with the seqid,
# length, and sequnce for each zotu.
outfile = open("pasteuria_zotu_metadata_rerun.csv", "w")

# Write column headers to the output file.
outfile.write("Zotu_ID,seqlen,seq\n")

# For each fasta record calculate the seq length and output to
# the csv file with the record id and sequence.
for title, seq in SimpleFastaParser(fasta_handle):
    seqlen = len(seq)
    outfile.write("%s,%s,%s\n" % (title, seqlen, seq))

outfile.close()
