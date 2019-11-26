#FASTQ filtering
#  --fastq_filter FILENAME     filter FASTQ file, output to FASTQ or FASTA file
# Options
#  --eeout                     include expected errors in FASTQ filter output
#  --fastaout FILENAME         FASTA output filename for nemsed sequences
#  --fastaout_discarded FNAME  FASTA filename for discarded sequences
#  --fastqout FILENAME         FASTQ output filename for nemsed sequences
#  --fastqout_discarded FNAME  FASTQ filename for discarded sequences
#  --fastq_ascii INT           FASTQ input quality score ASCII base char (33)
#  --fastq_maxee REAL          maximum expected error value for FASTQ filter
#  --fastq_maxee_rate REAL     maximum expected error rate for FASTQ filter
#  --fastq_maxns INT           maximum number of N's for FASTQ filter
#  --fastq_minlen INT          minimum length for FASTQ filter
#  --fastq_stripleft INT       bases on the left to delete for FASTQ filter
#  --fastq_trunclen INT        read length for FASTQ filter truncation
#  --fastq_truncqual INT       base quality value for FASTQ filter truncation
#  --relabel STRING            relabel filtered sequences with given prefix
#  --relabel_keep              keep the old label after the new when relabelling
#  --relabel_md5               relabel filtered sequences with md5 digest
#  --relabel_sha1              relabel filtered sequences with sha1 digest
#  --sizeout                   include abundance information when relabelling
#  --xsize                     strip abundance information in output

cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/QC-trimming-and-assembly

# Filter the concatenated reads to exclude those
# with an expected error > 1 using vsearch
# see https://doi.org/10.1093/bioinformatics/btv401
# for an explanation of minimum expected error filtering)

# Where possible, we have used vsearch
# because the underlying code is publicly available
vsearch --fastq_filter all-nem-amplicons.fastq \
--fastqout all-nem-amp-maxeefiltered.fastq \
--fastaout all-nem-amp-maxeefiltered.fasta \
--fastqout_discarded nem-failedexp-err-filter.fastq \
--fastq_maxee 1 --fastq_maxns 0

wait

vsearch --fastq_filter all-pas-amplicons.fastq \
--fastqout all-pas-amp-maxeefiltered.fastq \
--fastaout all-pas-amp-maxeefiltered.fasta \
--fastqout_discarded pas-failedexp-err-filter.fastq \
--fastq_maxee 1 --fastq_maxns 0
