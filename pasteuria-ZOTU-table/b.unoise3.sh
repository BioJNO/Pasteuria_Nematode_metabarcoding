cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table

# Instead of clustering amplicons based on sequence similarity we are going to
# try and capture all exact sequence variants (real biological sequences).
# See here for justification ("https://doi.org/10.1038/ismej.2017.119").
#
# To do this we're using unoise3 which is part of usearch ("https://doi.org/10.1101/081257").
#
# Briefly, this takes two sequences (A+B) and compares the number of positional differences
# between them and the abundance skew (abundance of A/abundance of B).
#
# This formula is used: B(d)=1/2^ad+1, where "d" is the number of positional
# differences and "a" is an alpha value set by the user
#
# If the abundance skew is less than B(d) then sequence A is judged to be the
# result of a sequencing/PCR error of exact sequence variant B; sequence A's
# abundance is added to sequence B's to make a Zero-radius OTU (ZOTU) with sequence=B.
# If not then both are judged to be real biological sequnces and ZOTUs A and B are created.
#
# Increasing the alpha value increases the value of B(d) and thus the chances of sequence
# A being assigned as an error.
#
# Here we use a low alpha value for as high a resolution as possible.
#
# For example two sequences, A and B, have an abundance skew of 1/5 and a
# levenshtein distance of 1. With a=2, B(1)= 1/8 so sequence A is NOT assigned
# as a mis-read of sequence B but is treated as a real biological sequence
# that may otherwise have been folded into a cluster with sequnce B. Thus
# zotus allow greater resolution of metabarcode sequences.

usearch -unoise3 ../QC-trimming-and-assembly/global-pas-derep.fasta \
-unoise_alpha 1 \
-tabbedout pasteuria_unoise3.txt \
-zotus pasteuria_alpha1_zotus.fasta

# The version of usearch we are using does not like ZOTU instead of OTU as
# the SeqID when it comes time to build the table so we will rename them.
sed -i -e 's/Zotu/Otu/g' pasteuria_alpha1_zotus.fasta

# Next we'll gather some ZOTU metadata into a tabular format
# which we will use later on to filter the final ZOTU table.
python b.extract_zotu_metadata.py