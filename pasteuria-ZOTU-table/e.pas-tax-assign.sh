cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table

# Now that we have our ZOTU table with the abundance of each ZOTU in each
# sample we're going to add metadata. We're interested in how each ZOTU
# compares to sequences related to organisms we already know something about.
#
# To do this comparison first we need a trustworthy database of sequences
# linked to taxonomic IDs. Such a database does not exist for Pasteuria spp.
# so we have created our own curated reference set from cloned sequences and
# from the publicly available sequences in the Silva-108 format:
# pas-refs.fasta, pas-taxmap.txt.
#
# Now we have our ZOTUs and trimmed reference database we can assign taxonomy
# to them. We use uclust and assign_taxonomy.py to do this which is part of
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
    assign_taxonomy.py -i pasteuria_alpha1_zotus.fasta \
                   -r pas-refs.fasta \
                   -t pas-taxmap.txt \
                   -o p$i-tax_match \
                   --similarity $i
    echo "done"

done

# We should now have a folder for each identity threshold with taxonomic
# assignments at that level. Each file within the foleders will have the
# same name "pasteuria_alpha1_zotus_tax_assignments.txt". We want to give
# them unique names 

mkdir combined_taxonomy 
# Using another bash loop to give each file a unique name and put them in a
# single folder
echo "renaming assigned taxonomy files"
for i in 1.0 0.995 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9
do 
    cp p$i-tax_match/pasteuria_alpha1_zotus_tax_assignments.txt combined_taxonomy/$i-tax-assignment.txt

done

cd combined_taxonomy/

# Make a list of tax assignment handles, order them from best highest
# stringency to lowest (this is important!)
ls *-tax-assignment.txt | sort -n -r > tax_assignment_handles.txt

# Run the combine_tax_assignments script to genetate a single combined
# best match taxonomic assignment for each ZOTU.
python ../e.combine_tax_assignments.py