cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/pasteuria-ZOTU-table

# All analysis from now on will use R (in R studio).
# First we'll make a folder for our project.
mkdir pasteuria_ZOTU_data_R

# Then we need the three tables we've generated.
cp pasteuria_alpha1_zotutab.txt pasteuria_ZOTU_data_R/
cp pasteuria_zotu_metadata.csv pasteuria_ZOTU_data_R/
cp combined_taxonomy/pasteuria_taxonomy_combined.txt pasteuria_ZOTU_data_R/

# Merging these tables will require a column with a unique header which is
# common to all tables.
#
# Both of the python generated tables have a Zotu_ID column but we need to
# rename the zotutab ID column
cd pasteuria_ZOTU_data_R
sed -i -e 's/#OTU ID/Zotu_ID/g' pasteuria_alpha1_zotutab.txt

# The IDs are still Otu1 etc but that doesn't matter because we have new,
# more informative, names to give them from the combined taxonomy table.

# I prefer to work with csv files in R so we convert the tsvs.
sed -i -e 's/\t/,/g' pasteuria_alpha1_zotutab.txt
sed -i -e 's/\t/,/g' pasteuria_taxonomy_combined.txt

# Then rename them.
mv pasteuria_alpha1_zotutab.txt pasteuria_alpha1_zotutab.csv
mv pasteuria_taxonomy_combined.txt pasteuria_taxonomy_combined.csv