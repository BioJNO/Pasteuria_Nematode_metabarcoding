cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/

# Make a new directory to work with pasteuria amplicon data and change into it
mkdir pasteuria-ZOTU-table
cd pasteuria-ZOTU-table

# Copy the Vsearch filtered assembled nematode reads to this directory (.)
cp ../QC-trimming-and-assembly/all-pas-amp-maxeefiltered.fastq .
# Wait for this to finish
wait

# Make a directory to store the amplicons sorted into samples
mkdir sample_fastq

# Invoke the sorting script 
python A.sort-pas-allowing-primer-mismatch.py
wait

# Move all the sample sorted fastq files into the directory you created.
mv pas*.fastq sample_fastq/