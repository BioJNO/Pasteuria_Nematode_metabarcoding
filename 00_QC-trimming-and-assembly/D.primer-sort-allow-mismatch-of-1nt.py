# Import from Biopython FastqGeneralIterator,
# which is a faster way to parse fastq files.
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Import regex which allows regular expression matching with mismatch.
import regex

# Write an output file to store all sorted assembled reads.
outfile1 = open("Pas-amplicons-containing-FandR-correct-orientation.fastq", "w")

# Sometimes about half the reads are reverse oriented,
# so we'll look for those too.
outfile2 = open("Pas-amplicons-containing-FandR-reverse-orientation.fastq", "w")
outfile3 = open("Nem-amplicons-containing-FandR-correct-orientation.fastq", "w")
outfile4 = open("Nem-amplicons-containing-FandR-reverse-orientation.fastq", "w")
outfile5 = open("Cyt-amplicons-containing-fandr-correct-orientation.fastq", "w")
outfile6 = open("Cyt-amplicons-containing-fandr-reverse-orientation.fastq", "w")
outfile7 = open("no-primer-pair-found.fastq", "w")

# Compile regular expressions (string) for each primer
# in each orientation, allowing a maximum mismatch of 1nt ({s<=1}).

PasF = regex.compile("(CAGCATCTTTGTGCCGAAGG){s<=1}")
PasR_rc = regex.compile("(TTGGAGAGACAGCCGGCG){s<=1}")

PasR = regex.compile("(CGCCGGCTGTCTCTCCAA){s<=1}")
PasF_rc = regex.compile("(CCTTCGGCACAAAGATGCTG){s<=1}")

NemF = regex.compile("(GGTGGTGCATGGCCGTTCTTAGTT){s<=1}")
NemR_rc = regex.compile("(ATTACGTCCCTGCCCTTTGTA){s<=1}")

NemR = regex.compile("(TACAAAGGGCAGGGACGTAAT){s<=1}")
NemF_rc = regex.compile("(AACTAAGAACGGCCATGCACCACC){s<=1}")

# These are primers used in a seperate study
# which shared the same sequncing run
CytF = regex.compile("(CTTGAAGACCTTCTGTAAAAATG){s<=1}")
CytR_rc = regex.compile("(CTCTTAAGACGGTAGCTCG){s<=1}")

CytR = regex.compile("(CGAGCTACCGTCTTAAGAG){s<=1}")
CytF_rc = regex.compile("(CATTTTTACAGAAGGTCTTCAAG){s<=1}")

# Open and store the input file handle as a variable.
input_handle = open("pear_assembled-default-q22-minlen150.assembled.fastq")

# For each fastq record in the assembled read file...
for title, seq, qual in FastqGeneralIterator(input_handle):
    # search for each primer pair in the assembled reads in forward...
    if PasF.search(seq) and PasR_rc.search(seq):
        outfile1.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    # and, if not, reverse orientation.
    elif PasR.search(seq) and PasF_rc.search(seq):
        outfile2.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    elif NemF.search(seq) and NemR_rc.search(seq):
        outfile3.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    elif NemR.search(seq) and NemF_rc.search(seq):
        outfile4.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    elif CytF.search(seq) and CytR_rc.search(seq):
        outfile5.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    elif CytR.search(seq) and CytF_rc.search(seq):
        outfile6.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    # If no primer pairs are found in any orientation
    # spit them out into this file so that you can check
    # sorting is working as you want it to.
    else:
        outfile7.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
