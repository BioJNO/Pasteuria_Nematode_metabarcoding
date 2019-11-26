cd /home/jo42324/scratch/metagenetics/Pasteuria_Nematode_metabarcoding/nematode-ZOTU-table

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

usearch -unoise3 sample_fastq/global-nem-derep.fasta \
-unoise_alpha 1 \
-tabbedout minee2-nematode_unoise3.txt \
-zotus minee2-nematode_alpha1_zotus.fasta

# The version of usearch we are using does not like ZOTU instead of OTU as
# the SeqID when it comes time to build the table so we will rename them.
sed -i -e 's/Zotu/Otu/g' minee2-nematode_alpha1_zotus.fasta

# Next we'll gather some ZOTU metadata into a tabular format
# (we use this later on to filter the final ZOTU table).
python e.extract_zotu_metadata.py
