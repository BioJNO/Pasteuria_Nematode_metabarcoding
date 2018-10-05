
cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table/sample_fastq

# Make a list of sorted fastq files as lines in a text file.
ls > pas-fq-files.txt

# This list will include the file you used to create it so remove that with sed.
sed -i '/pas-fq-files.txt/d' pas-fq-files.txt 

# Invoke the renaming python script to give each sequence a sample number
# for usearch and join them all to a single file.
python ../c.rename_fastq_samples.py
wait 

mv all-pas-fq-w-usearch-sampleids.fastq ../
