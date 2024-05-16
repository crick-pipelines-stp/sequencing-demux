#!/bin/bash

# Requires pod5 tools to be installed
# pip install pod5

# Check if the correct number of arguments is passed
if [ $# -ne 2 ]; then
    echo "Usage: $0 <path-to-directory> <total-number-of-reads>"
    exit 1
fi

# Assign command line arguments to variables
directory=$1
total_reads=$2

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Directory does not exist: $directory"
    exit 1
fi

# Count the number of files in the directory
file_count=$(ls -1 "$directory"/*.pod5 | wc -l)

# Check if there are any files
if [ "$file_count" -eq 0 ]; then
    echo "No .pod5 files found in the directory."
    exit 1
fi

# Calculate the number of reads per file
reads_per_file=$((total_reads / file_count))
echo "Files in dir = $file_count"
echo "Read Ids per file = $reads_per_file"

# Create a directory for subset and filtered files
subset_dir="$directory/subsets"
mkdir -p "$subset_dir"

# Iterate over each file in the directory
for file in "$directory"/*.pod5
do
    echo "Processing $file..."

    # Generate a subset file name
    base_name=$(basename "$file" .pod5)
    subset_reads="$subset_dir/${base_name}_subset_readids.txt"
    subset_file="$subset_dir/${base_name}_subset_$reads_per_file.pod5"

    # Extract read IDs using pod5 inspect, shuffle them, and get the first N reads_per_file
    pod5 inspect reads "$file" | awk -F, '{print $1}' | gshuf | head -n "$reads_per_file" > "$subset_reads"

    # Filter the pod5 based on the read ids
    pod5 filter "$file" --output "$subset_file" --ids "$subset_reads"
done

echo "complete!"
