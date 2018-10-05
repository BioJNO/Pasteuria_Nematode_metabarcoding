# Import required modules.
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex
from itertools import product
from lru import LRU

# Store primers as variables.
NemF = "GGTGGTGCATGGCCGTTCTTAGTT"
NemR = "TACAAAGGGCAGGGACGTAAT"
NemRseq = Seq(NemR, generic_dna)
NemR_rc = str(NemRseq.reverse_complement())

# Compile primer regular expressions allowing a positional mismatch
# (substitution) of one or less.
NemFs = regex.compile("(GGTGGTGCATGGCCGTTCTTAGTT){s<=1}")
NemR_rcs = regex.compile("(ATTACGTCCCTGCCCTTTGTA){s<=1}")

# Store list of barcode sequences as variable.
taglist = ["AAGGTC", "ACCTCA", "ACGTGT", "ACTCTG", "AGCATG", "AGTCCA",
           "CAACTC", "CAAGCA", "CACAGT", "CAGGAT", "CAGTTG", "CCATAC",
           "CCTGTA", "CGATCT", "CGTAGA", "CTCACA", "CTGAAC", "CTTGCT",
           "GAAGTG", "GACTTC", "GAGCTA", "GATGGT", "GCAGAA", "GCTTGA",
           "GGAACA", "GGTATC", "GTCGTA", "GTGATG", "GTTCAC", "TCCAGA",
           "TCGTTC", "TGGACT"]

# Split forward tags into 4 sets of 8 corresponding to primer dilution plates.
taglist1 = taglist[0:8]
taglist2 = taglist[8:16]
taglist3 = taglist[16:24]
taglist4 = taglist[24:32]

# Create an empty list to store reverse complement of barcodes.
rctaglist = []

# Using Biopython generate a list of reverse complement barcodes.
for tag in taglist:
    # Convert string to seq object.
    dna_tag = Seq(tag, generic_dna)
    # Genereate reverse complement of seq object.
    rctagseq = dna_tag.reverse_complement()
    # Convert reverse complement seq objects back to string.
    rctagstring = str(rctagseq)
    # Append rctags as string to list.
    rctaglist.append(rctagstring)

# Split reverse comp tags into 4 sets of 8 corresponding to dilution plates.
rctaglist1 = rctaglist[0:8]
rctaglist2 = rctaglist[8:16]
rctaglist3 = rctaglist[16:24]
rctaglist4 = rctaglist[24:32]

# Generate a list of all forward and reverse barcode combinations.
frcombos1 = list(product(rctaglist1, taglist1))
frcombos2 = list(product(rctaglist2, taglist1))
frcombos3 = list(product(rctaglist3, taglist1))
frcombos4 = list(product(rctaglist4, taglist1))

frcombos5 = list(product(rctaglist1, taglist2))
frcombos6 = list(product(rctaglist2, taglist2))
frcombos7 = list(product(rctaglist3, taglist2))
frcombos8 = list(product(rctaglist4, taglist2))

frcombos9 = list(product(rctaglist1, taglist3))
frcombos10 = list(product(rctaglist2, taglist3))
frcombos11 = list(product(rctaglist3, taglist3))
frcombos12 = list(product(rctaglist4, taglist3))

frcombos13 = list(product(rctaglist1, taglist4))
frcombos14 = list(product(rctaglist2, taglist4))
frcombos15 = list(product(rctaglist3, taglist4))
frcombos16 = list(product(rctaglist4, taglist4))

fr_combo_final = (frcombos1 + frcombos2 + frcombos3 + frcombos4 +
                  frcombos5 + frcombos6 + frcombos7 + frcombos8 +
                  frcombos9 + frcombos10 + frcombos11 + frcombos12 +
                  frcombos13 + frcombos14 + frcombos15 + frcombos16)

# Define sample file handle
input_handle = open("all-nem-amp-maxeefiltered.fastq")

# We need a file to store assembled reads from each sample.

# There is a system imposed limit on the number of files we can have open at
# any one time so here we use a last recently used cache.
#
# At most 999 files open at one time from the LRU cache.
cache = LRU(999)

# Create a dictionary to store the output filenames.
filenames = {}
# Create an output fq file for each sample which is numbered, and contains the
# corresponding barcode pair
for x, (rbar, fbar) in enumerate(fr_combo_final):
    f = "nem-%05i-%s-%s.fastq" % (x, fbar, rbar)
    filenames[(rbar, fbar)] = f
    outfile = open(f, "w")
    outfile.close()

# Create an output file to store amplicons which have primer sequences but
# none of the expected barcode pairs
oddities = open('nem-oddities.txt', 'w')

# Create an empty numeric to store a tally of fragmented barcodes
fragmentary = 0
# Create an empty dictionary for expected barcode pair tallies
fr_tallies = dict()

# Using FastqGeneralIterator (because it is a faster way to parse fq files)
# search each sequence for the 6nt seq before and after primers.
for title, seq, qual in FastqGeneralIterator(input_handle):
    seqlen = len(seq)
    fprimersear = NemFs.search(seq)
    rprimersear = NemR_rcs.search(seq)
    fstart = fprimersear.start()
    rstart = rprimersear.start()
    fend = fprimersear.end()
    rend = rprimersear.end()
    frame = seq[fstart:rend]
    fbar = seq[fstart - 6: fstart]
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
        if (fbar, rbar) in fr_tallies:
            # The a+=b trick is short for a=a+b
            # The notation comes from the C langauge.
            fr_tallies[(fbar, rbar)] += 1
        else:
            # First time to see it, count it
            fr_tallies[(fbar, rbar)] = 1
        # Do we already have an output file ready for this pair?
        if (rbar, fbar) in filenames and (rbar, fbar) in cache:
            # Yes, write the frame to file
            print("Wanted  %s %s" % (fbar, rbar))
            cache[(rbar, fbar)].write("@%s\n%s\n+\n%s\n" % (title, noprimframe, noprimqual))
        elif (rbar, fbar) in filenames:
            # No, open the file from the dictionary of handles
            name = filenames[(rbar, fbar)]
            cache[(rbar, fbar)] = open(name, "a")
            # Then write the frame to file
            cache[(rbar, fbar)].write("@%s\n%s\n+\n%s\n" % (title, noprimframe, noprimqual))
        else:
            # We didn't put this barcode pair here!
            # Write it to the oddities file.
            print("Unexpected %s %s" % (fbar, rbar))
            oddities.write("%s\t%s\t%s\t%s\n" % (fbar, rbar, title, seq))

# Close open files.
oddities.close()
input_handle.close()

# Print the tallies for each expected barcode pair.
print("Observed barcode pairs, tally count, wanted or not?")
for (fbar, rbar) in fr_tallies:
    print("%s %s count %i %r" % (fbar, rbar, fr_tallies[(fbar, rbar)], (fbar, rbar) in fr_combo_final))

# Print out the full length unexpected barcode and fragment tallies.
print("In total %i full length barcode pairs, and %i fragments" % (sum(fr_tallies.values()), fragmentary))
