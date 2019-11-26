cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/QC-trimming-and-assembly

# Invoke the correct orientation python script.
python E.correct-orientation.py

# Wait for it to finish.
wait 

# Concatenate flipped reads and reads which assembled in the expected orientation.
cat pas-amplicons-orientation-corrected.fastq \
Pas-amplicons-containing-FandR-correct-orientation.fastq \
> all-pas-amplicons.fastq

wait

# For both primer sets.
cat nem-amplicons-orientation-corrected.fastq \
Nem-amplicons-containing-FandR-correct-orientation.fastq \
> all-nem-amplicons.fastq
