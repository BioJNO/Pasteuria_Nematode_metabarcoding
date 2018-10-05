cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/nematode-ZOTU-table

# Map fq files to zotus to generate a ZOTU table.
usearch -otutab all-nem-fq-w-usearch-sampleids.fastq \
-otus nematode_alpha1_zotus.fasta \
-otutabout nematode_alpha1_zotutab.txt \
-mapout nematode_alpha1_zotu_map.txt \
-notmatched nematode_alpha1_unmapped.fasta
