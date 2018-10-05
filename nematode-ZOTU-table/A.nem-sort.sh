cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/

# Make a new directory to work with nematode amplicon data and change into it
mkdir nematode-ZOTU-table
cd nematode-ZOTU-table

# Copy the Vsearch filtered assembled nematode reads to this directory (.)
cp ../QC-trimming-and-assembly/all-nem-amp-maxeefiltered.fastq .
# Wait for this to finish
wait

# Make a directory to store the amplicons sorted into samples
mkdir sample_fastq

# Invoke the sorting script 
python A.sort-nem-allowing-primer-mismatch.py
wait

# Move all the sample sorted fastq files into the directory you created.
mv nem*.fastq sample_fastq/
