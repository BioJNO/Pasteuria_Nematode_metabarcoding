cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/nematode-ZOTU-table

# Now that we have our ZOTU table with the abundance of each ZOTU in each
# sample we're going to add metadata. We're interested in how each ZOTU
# compares to sequences related to organisms we already know something about.
#
# To do this comparison first we need a trustworthy database of sequences
# linked to taxonomic IDs. For our nematode dataset we're using the Silva-132
# Nr99 database ("https://www.arb-silva.de/"), a curated database which
# covers our amplification target (18S SSU).
# 
# We've trimmed this to the amplified region between primers and removed
# uncultred environmental samples as uninformative
#
# The resultant trimmed Silva-132 database is available on Figshare:
#
# https://doi.org/10.6084/m9.figshare.9897374.v2 
#
# We use uclust and assign_taxonomy.py to do taxonomic assignments via
# qiime.
#
# Get qiime using (conda install -c bioconda qiime)
# We need to activate the qiime environment first
source activate qiime1

# Now we want to assign taxonomy iteratively so that we can get the best
# possible match for each sequence given that sequences are not clustered
# and we are considering them real, biologically relevant, sequence variants.
# We'll use a bash loop to do this which outputs a folder for each identity
# threshold.
# For 100-90% identity...
for i in 1.0 0.995 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9
do 
    echo "running similarity threshold $i taxa mapping"
    assign_taxonomy.py -i nematode_alpha1_zotus.fasta \
                   -r /home/jo42324/scratch/silva/reduced_ref_sequences.fasta \
                   -t /home/jo42324/scratch/silva/reduced_ref_taxonomy.txt \
                   -o p$i-tax_match \
                   --similarity $i
    echo "done"

done

# We should now have a folder for each identity threshold with taxonomic
# assignments at that level. Each file within the folders will have the
# same name "nematode_alpha1_zotus_tax_assignments.txt". We want to give
# them unique names 

mkdir combined_taxonomy 
# Using another bash loop to give each file a unique name and put them in a
# single folder
echo "renaming assigned taxonomy files"
for i in 1.0 0.995 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9
do 
    cp p$i-tax_match/nematode_alpha1_zotus_tax_assignments.txt combined_taxonomy/$i-tax-assignment.txt

done

cd combined_taxonomy/

# Make a list of tax assignment handles, order them from best highest
# stringency to lowest (this is important!)
ls *-tax-assignment.txt | sort -n -r > tax_assignment_handles.txt

# Run the combine_tax_assignments script to genetate a single combined
# best match taxonomic assignment for each ZOTU.
python ../f.combine_tax_assignments.py
