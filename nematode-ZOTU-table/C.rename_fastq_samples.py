from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Open the list of fq file handles
handle_list = open("nem-fq-files.txt")
# Write an output file to hold all fq records
f_out = open("all-nem-fq-w-usearch-sampleids.fastq", "w")

# Set an empty numeric value to hold the sample number as we iterate through the sample fq files.
n = 0

# For each sample fq file in the list
for handle in handle_list:
    # Get the handle name and strip the newline character
    file_name = handle.rstrip("\n")
    # Open the sample fq file
    input_handle = open(file_name)
    # Update the sample number
    n = n+1
    for title, seq, qual in FastqGeneralIterator(input_handle):
        # Set a new name based on the old one which is acceptable to usearch
        newname = title + ";sample=" + str(n) + ";"
        # Write it to the combined output file
        f_out.write("@%s\n%s\n+\n%s\n" % (newname, seq, qual))

# Close the output file.
f_out.close()
