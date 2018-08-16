from Bio.SeqIO.FastaIO import SimpleFastaParser
import regex

# Store primer regular expressions as variables
NemF = regex.compile("(GGTGGTGCATGGCCGTTCTTAGTT){s<=1}")
# This is reverse compliment!
NemR_rc = regex.compile("(ATTACGTCCCTGCCCTTTGTA){s<=1}")

# Store the Silva 108 eukaryotic fasta set as a variable
# We've added plasmid sequneces we generated for the controls to this database and to the taxonomy.txt file.
input_handle = open("Silva_108_rep_set_Eukarya_only.fna")

# Open an output file for each condition test in the cutprimers function
outfile1 = open("Silva-108-primer-trim-strict-len.fasta", "w")
outfile2 = open("Silva-108-failed-minlen.fasta.fasta", "w")
outfile3 = open("Silva-108-failed-primer-maxlen-filter.fasta", "w")
outfile4 = open("Silva-108-no-primer-match.fasta", "w")


# We've set this up as a function to make it easier to fiddle with the parameters
def cutprimers(input_hadle, fprimer, rprimer, minlen, maxlen, outfile1, outfile2, outfile3, outfile4):
	# Parse the fasta reference file
	for title, seq in SimpleFastaParser(input_handle):
		# Look for the f and r primers in each sequence
		fprimersear = fprimer.search(seq)
		rprimersear = rprimer.search(seq)

		# If you've found both primers...
		if fprimersear > -1 and rprimersear > -1:
			# Get the start and end position of the primers
			fstart = fprimersear.start()
			rstart = rprimersear.start()
			fend = fprimersear.end()
			rend = rprimersear.end()
			# Get the frame, which will be the size of your PCR products without barcodes
			frame = seq[fstart:rend]
			# Get the frame between primers so that proportional variation between sequences is maximised
			noprimframe = seq[fend:rstart]
			# Calculte the thoretical length of the PCR product using these primers
			seqlen = len(frame)

			# if the size is correct...
			if seqlen >= minlen and seqlen <= maxlen:
				outfile1.write(">%s\n%s\n" % (title, noprimframe))

			# If the product is too short...
			elif seqlen < minlen:
				outfile2.write(">%s\n%s\n" % (title, noprimframe))

			# If the product is too long...
			elif seqlen > maxlen:
				outfile3.write(">%s\n%s\n" % (title, noprimframe))

		# If you dont find the primers...
		else:
			outfile4.write(">%s\n%s\n" % (title, seq))


# Call the function
cutprimers(input_handle, NemF, NemR_rc, 339, 420, outfile1, outfile2, outfile3, outfile4)

# Make sure all the output files are closed
outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()

# Now we need to trim the taxa mapping file that goes with the Silva Database

# Open the trimmed fasta database
Trimmed_silva_fa = open("Silva-108-primer-trim-strict-len.fasta", "r")
# Open the untrimmed taxa mapping file
Silva_tax_map = open("Silva_108_taxa_mapping.txt", "r")
# Write a file to store the trimmed taxa mapping data
Silva_tax_map_trimmed = open("Silva_108_taxa_mapping_trimmed.txt", "w")

# Create an empty dictionary of sequence names
name_dict = {}
not_matched = open("not_matched_tax.txt", "w")

# For each line in the taxa mapping file
for line in Silva_tax_map:
	# get the name and store it as a variable
	namestop = line.find("\t")
	name = line[:namestop]
	# get the taxonomic info and store it as a variable
	tax = line[namestop:]
	# store the name and the taxonomic assignment as a dictionary item
	f = name, tax
	name_dict[(name, tax)] = f

# For each sequence in the trimmed fasta file
for title, seq in SimpleFastaParser(Trimmed_silva_fa):
	# Iterate through the dictionary of ID:taxonomic_info
	for (name, tax) in name_dict:
		# if the sequnce title matches the name in the dictionary...
		if name == title:
			# write the sequence ID and tax record to the trimmed taxa file
			name = name.rstrip("\n")
			Silva_tax_map_trimmed.write("%s%s\n" % (name, tax))
		else:
			name = name.rstrip("\n")
			not_matched.write("%s%s\n" % (name, tax))
