cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/QC-trimming-and-assembly

# Assemble trimmed PE illumina reads with PEAR.
/home/jo42324/Download/PEAR-master/bin/pear \
-f trim-q22-minlen150/10623_metagen_trimmed_f.fq \
-r trim-q22-minlen150/10623_metagen_trimmed_r.fq \
-j 16 -o pear_assembled-default-q22-minlen150
