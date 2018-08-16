# Last we collect some sequence metadata which we'll use in data filtering

from Bio.SeqIO.FastaIO import SimpleFastaParser

fasta_handle = open("nematode_alpha1_zotus.fasta")

outfile = open("nematode_zotu_metadata.csv", "w")
# Write column headers to your output file
# We're going to make a csv table with the seqid, length, and sequnce for each zotu
outfile.write("Zotu_ID,seqlen,seq\n")

for title, seq in SimpleFastaParser(fasta_handle):
    seqlen = len(seq)
    outfile.write("%s,%s,%s\n" % (title, seqlen, seq))

outfile.close()
