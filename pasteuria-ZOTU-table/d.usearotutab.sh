cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table

# Map fq files to zotus to generate a ZOTU table.
usearch -otutab all-pas-fq-w-usearch-sampleids.fastq \
-otus pasteuria_alpha1_zotus.fasta \
-otutabout pasteuria_alpha1_zotutab.txt \
-mapout pasteuria_alpha1_zotu_map.txt \
-notmatched pasteuria_alpha1_unmapped.fasta
