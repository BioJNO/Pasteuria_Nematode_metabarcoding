
cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/nematode-ZOTU-table/sample_fastq

# Make a list of sorted fastq files as lines in a text file.
ls > nem-fq-files.txt

# This list will include the file you used to create it so remove that with sed.
sed -i '/nem-fq-files.txt/d' nem-fq-files.txt 

# Move the list back up one level 
mv list.txt ../

# Follow it 
cd ../

# Invoke the renaming python script to give each sequence a sample number
# for usearch and join them all to a single file.
python C.rename_fastq_samples.py

