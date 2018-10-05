#!/usr/bin/env python

# Import required modules
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import regex
from itertools import product
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from lru import LRU

# Store primers as string variables
PasF = "CAGCATCTTTGTGCCGAAGG"
PasR = "CGCCGGCTGTCTCTCCAA"

# Convert reverse primer sequence string to seq object,
# then convert to reverse complement and store this as string
PasRseq = Seq(PasR, generic_dna)
PasR_rc = str(PasRseq.reverse_complement())

# Compile primers as regular expressions allowing
# a positional mismatch (substitution) of one or less
PasFs = regex.compile("(CAGCATCTTTGTGCCGAAGG){s<=1}")
PasR_rcs = regex.compile("(TTGGAGAGACAGCCGGCG){s<=1}")

# Store list of barcode sequences as variable

taglist = ["AAGGTC", "ACCTCA", "ACGTGT", "ACTCTG", "AGCATG", "AGTCCA",
           "CAACTC", "CAAGCA", "CACAGT", "CAGGAT", "CAGTTG", "CCATAC",
           "CCTGTA", "CGATCT", "CGTAGA", "CTCACA", "CTGAAC", "CTTGCT",
           "GAAGTG", "GACTTC", "GAGCTA", "GATGGT", "GCAGAA", "GCTTGA",
           "GGAACA", "GGTATC", "GTCGTA", "GTGATG", "GTTCAC", "TCCAGA",
           "TCGTTC", "TGGACT"]

# Split forward tags into 4 sets of 8 corresponding to primer dilution plates
taglist1 = taglist[0:8]
taglist2 = taglist[8:16]
taglist3 = taglist[16:24]
taglist4 = taglist[24:32]

# Create an empty list to store reverse complement of barcodes
rctaglist = []

# Using Biopython generate a list of reverse complement barcodes
for tag in taglist:
    # Convert string to seq object
    dna_tag = Seq(tag, generic_dna)
    # Use .reverse_complement function to genereate
    # reverse complement of seq object
    rctagseq = dna_tag.reverse_complement()
    # Convert reverse complement seq objects back to string
    rctagstring = str(rctagseq)
    # Append rctags as string to list
    rctaglist.append(rctagstring)

# Split reverse complement barcodes into
# 4 sets of 8 corresponding to dilution plates
rctaglist1 = rctaglist[0:8]
rctaglist2 = rctaglist[8:16]
rctaglist3 = rctaglist[16:24]
rctaglist4 = rctaglist[24:32]

# Generate a list of all forward and reverse barcode combinations
frcombos1 = list(product(taglist1, rctaglist1))
frcombos2 = list(product(taglist2, rctaglist1))
frcombos3 = list(product(taglist3, rctaglist1))
frcombos4 = list(product(taglist4, rctaglist1))

frcombos5 = list(product(taglist1, rctaglist2))
frcombos6 = list(product(taglist2, rctaglist2))
frcombos7 = list(product(taglist3, rctaglist2))
frcombos8 = list(product(taglist4, rctaglist2))

frcombos9 = list(product(taglist1, rctaglist3))
frcombos10 = list(product(taglist2, rctaglist3))
frcombos11 = list(product(taglist3, rctaglist3))
frcombos12 = list(product(taglist4, rctaglist3))

frcombos13 = list(product(taglist1, rctaglist4))
frcombos14 = list(product(taglist2, rctaglist4))
frcombos15 = list(product(taglist3, rctaglist4))
frcombos16 = list(product(taglist4, rctaglist4))

fr_combo_final = (frcombos1 + frcombos2 + frcombos3 + frcombos4 +
                  frcombos5 + frcombos6 + frcombos7 + frcombos8 +
                  frcombos9 + frcombos10 + frcombos11 + frcombos12 +
                  frcombos13 + frcombos14 + frcombos15 + frcombos16)

# Define sample file handle
input_handle = open("all-pas-amp-maxeefiltered.fastq")

# Open all the output handles at the start, and keep them in a dictionary
cache = LRU(999)

filenames = {}
for x, (fbar, rbar) in enumerate(fr_combo_final):
    f = "pas-%05i-%s-%s.fastq" % (x, fbar, rbar)
    filenames[(fbar, rbar)] = f
    outfile = open(f, "w")
    outfile.close()

oddities = open('pas-oddities.txt', 'w')

fragmentary = 0
fr_tallies = dict()

# Using SimpleFastaParser search each sequence
# for the 6nt seq before and after primers
for title, seq, qual in FastqGeneralIterator(input_handle):
    seqlen = len(seq)
    fprimersear = PasFs.search(seq)
    rprimersear = PasR_rcs.search(seq)
    fstart = fprimersear.start()
    rstart = rprimersear.start()
    fend = fprimersear.end()
    rend = rprimersear.end()
    frame = seq[fstart:rend]
    fbar = seq[fstart-6:fstart]
    rbar = seq[rend:rend + 6]
    # Store the amplified region between primers
    noprimframe = seq[fend:rstart]
    noprimqual = qual[fend:rstart]
    # Hoping to find pair of 6bp known barcodes.
    #
    # Checking against the expected pairs via the dictionary ensures
    # will only write this out once, without needing a for loop.
    #
    # Might have only partial sequences, e.g. 'CTGA' and 'GG'
    # Might have pcr or sequencing erors returning barcodes not on our list
    if len(fbar) != 6 or len(rbar) != 6:
        print("Ignoring %s %s" % (fbar, rbar))
        fragmentary += 1
    else:
        # Right length, first count the barcode pair
        # TODO: Use try/except with KeyError
        # TODO: Replace with a default dictionary?
        # Alternative style:
        # fr_tallies[(fbar, rbar)] = fr_tallies.get((fbar, rbar), 0) + 1
        if (fbar, rbar) in fr_tallies:
            # The a+=b trick is short for a=a+b
            # The notation comes from the C langauge.
            fr_tallies[(fbar, rbar)] += 1
        else:
            # First time to see it, count it
            fr_tallies[(fbar, rbar)] = 1
        if (fbar, rbar) in filenames and (fbar, rbar) in cache:
            print("Wanted  %s %s" % (fbar, rbar))
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        elif (fbar, rbar) in filenames:
            name = filenames[(fbar, rbar)]
            cache[(fbar, rbar)] = open(name, "a")
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        else:
            print("Unexpected %s %s" % (fbar, rbar))
            oddities.write("%s\t%s\t%s\t%s\n" % (fbar, rbar, title, seq))
        # Do we already have an output file ready for this pair?
        # There is a limit to how many files we can open at once...
        # if (fbar, rbar) not in out_handles:
        #    out_handles[(fbar, rbar)] =
        #    open("Pas-" + fbar + "-" + rbar + ".fasta", "w")
        # out_handles[(fbar, rbar)].write(">%s\n%s\n" % (title, seq))

#    Short cut for testing:
#    if sum(fr_tallies.values()) > 1000:
#        break

oddities.close()

input_handle.close()

# Print the tallies for each expected barcode pair
print("Observed barcode pairs, tally count, wanted or not?")
for (fbar, rbar) in fr_tallies:
    print("%s %s count %i %r" %
          (fbar, rbar, fr_tallies[(fbar, rbar)],
           (fbar, rbar) in fr_combo_final))

print("In total %i full length barcode pairs, and %i fragments" %
      (sum(fr_tallies.values()), fragmentary))
