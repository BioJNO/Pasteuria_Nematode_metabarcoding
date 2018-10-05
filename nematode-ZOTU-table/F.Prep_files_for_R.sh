cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/nematode-ZOTU-table

# All analysis from now on will use R (in R studio).
# First we'll make a folder for our project.
mkdir nematode_ZOTU_data_R

# Then we need the three tables we've generated.
cp nematode_alpha1_zotutab.txt nematode_ZOTU_data_R/
cp nematode_zotu_metadata.csv nematode_ZOTU_data_R/
cp combined_taxonomy/nematode_taxonomy_combined.txt nematode_ZOTU_data_R/

# Merging these tables will require a column with a unique header which is
# common to all tables.
#
# Both of the python generated tables have a Zotu_ID column but we need to
# rename the zotutab ID column
sed -i -e 's/#OTU ID/Zotu_ID/g' nematode_alpha1_zotutab.txt

# The IDs are still Otu1 etc but that doesn't matter because we have new,
# more informative, names to give them from the combined taxonomy table.

# I prefer to work with csv files in R so we convert the tsvs.
sed -i -e 's/\t/,/g' nematode_alpha1_zotutab.txt
sed -i -e 's/\t/,/g' nematode_taxonomy_combined.txt

# Then rename them.
mv nematode_alpha1_zotutab.txt nematode_alpha1_zotutab.csv
mv nematode_taxonomy_combined.txt nematode_taxonomy_combined.csv