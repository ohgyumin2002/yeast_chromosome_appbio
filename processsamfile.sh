#!/bin/bash

# Define the list of SAM files that we will process
FILES=(
    "sy14_r1_1.sam" "sy14_r1_2.sam" "sy14_r2_1.sam" "sy14_r2_2.sam"
    "sy14_r3_1.sam" "sy14_r3_2.sam" "wt_r1_1.sam" "wt_r1_2.sam"
    "wt_r2_1.sam" "wt_r2_2.sam" "wt_r3_1.sam" "wt_r3_2.sam"
)

# Path to annotation file (modify this path to your actual annotation file)
ANNOTATION="annotation_sam_file_processing.gtf"

# Create output directory for processed files
mkdir -p processed

# Step 1: Process each SAM file
for FILE in "${FILES[@]}"; do
    # Extract the base name from the file name
    BASENAME=$(basename "$FILE" .sam)

    echo "Processing $FILE ..."

    # Convert SAM to BAM
    samtools view -Sb "$FILE" > "processed/$BASENAME.bam"

    # Sort BAM
    samtools sort "processed/$BASENAME.bam" -o "processed/$BASENAME.sorted.bam"

    # Index BAM
    samtools index "processed/$BASENAME.sorted.bam"

    echo "$FILE converted, sorted, and indexed."
done

# Step 2: Count reads using featureCounts
SORTED_BAMS=$(find processed -name "*.sorted.bam" | tr '\n' ' ')

# Output count file
COUNTS_OUTPUT="gene_counts.txt"

featureCounts -a "$ANNOTATION" -o "$COUNTS_OUTPUT" $SORTED_BAMS

echo "Read counting complete. Results saved to $COUNTS_OUTPUT."

