#!/usr/bin/env bash

# Function to display usage information
usage() {
    echo "Usage: $0 <fasta_file> <tsv_file> <output_file>"
    echo "Options:"
    echo "  -h    Display this help message"
}

# Check for -h flag
if [ "$1" == "-h" ]; then
    usage
    exit 0
fi

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

# Define file paths from arguments
fasta_file="$1"
tsv_file="$2"
output_file="$3"

# Read protein names from the TSV file into an array
mapfile -t protein_names < <(cut -f1 "$tsv_file")

# Initialize a temporary file to store protein names
temp_file=$(mktemp)

# Write protein names to the temporary file
for protein in "${protein_names[@]}"; do
    echo "$protein" >>"$temp_file"
done

# Use seqkit to extract the sequences
seqkit grep -f "$temp_file" "$fasta_file" >"$output_file"

# Remove the temporary file
rm "$temp_file"

echo "Protein extraction completed. Extracted proteins are saved in $output_file."
