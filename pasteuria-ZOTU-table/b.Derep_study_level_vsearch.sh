# Vsearch dereplication

# Now that we have sorted and filtered our assembled
# reads (amplicons) we need to take all the amplicons
# binned (sorted into files for each primer pair)
# in the previous steps and de-replicate them at the
# study level (collapse all identical amplicons into
# a single record annotated with the total number of
# times it occurs).

# change into the directory where your binned and
# filtered fastq files are.
cd /home/jo42324/scratch/metagenetics/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table/sample_fastq

# Make join sorted sample files into a global reference set
cat pas-*.fastq > global-pas.fq

vsearch --derep_fulllength global-pas.fq \
--output global-pas-derep.fasta \
--sizeout
