# Now that we've sorted the amplicons into primer sets
# we need to flip the amplicons which assemble in the
# reverse orientation back into the expected orientation.


# Import the required modules from Biopython.
from Bio import SeqIO

# Open output files to store the flipped amplicons.
outfile1 = open("pas-amplicons-orientation-corrected.fastq", "w")
outfile2 = open("nem-amplicons-orientation-corrected.fastq", "w")
outfile3 = open("cyt-amplicons-orientation-corrected.fastq", "w")

# For each fastq record in the reverse
# orientation output from the previous step...
for seq_record in SeqIO.parse("Pas-amplicons-containing-FandR-reverse-orientation.fastq", "fastq"):
    # ...flip it. Keeping the name and appending "(RC)"
    # to the title so that later on we can see it was flipped.
    rcseq = seq_record.reverse_complement(id=True, description=seq_record.description + " (RC)")
    # Then write it to the output file.
    SeqIO.write(rcseq, outfile1, "fastq")

for seq_record in SeqIO.parse("Nem-amplicons-containing-FandR-reverse-orientation.fastq", "fastq"):
    rcseq = seq_record.reverse_complement(id=True, description=seq_record.description + " (RC)")
    SeqIO.write(rcseq, outfile2, "fastq")

for seq_record in SeqIO.parse("Cyt-amplicons-containing-fandr-reverse-orientation.fastq", "fastq"):
    rcseq = seq_record.reverse_complement(id=True, description=seq_record.description + " (RC)")
    SeqIO.write(rcseq, outfile3, "fastq")

# Close all the output files when you're done.
outfile1.close()
outfile2.close()
outfile3.close()
