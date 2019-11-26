cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/QC-trimming-and-assembly

# Run trimmomatic on raw PE Illumina reads
# setting a quality score cutoff of 22
# and a minimum length of 150
java -jar /home/jo42324/Download/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 4 \
-phred33 /home/jo42324/metabarcode/data/raw_data/Raw_illumina/10623_metagen_f_reads.fastq.gz \
/home/jo42324/metabarcode/data/raw_data/Raw_illumina/10623_metagen_r_reads.fastq.gz \
trim-q22-minlen150/10623_metagen_trimmed_f.fq \
trim-q22-minlen150/unpaired_metagen_f.fq \
trim-q22-minlen150/10623_metagen_trimmed_r.fq \
trim-q22-minlen150/unpaired_metagen_r.fq \
LEADING:22 TRAILING:22 SLIDINGWINDOW:10:22 MINLEN:150
