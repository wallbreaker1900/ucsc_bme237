#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH --partition=256x44
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Define the list of BAM files
BAM_FILES=(
    "SRR13890278_HISAT2_sorted.bam"
    "SRR13890277_HISAT2_sorted.bam"
    "SRR13890275_HISAT2_sorted.bam"
    "SRR13890274_HISAT2_sorted.bam"
    "SRR13890272_HISAT2_sorted.bam"
    "SRR13890271_HISAT2_sorted.bam"
    "SRR13890279_HISAT2_sorted.bam"
    "SRR13890276_HISAT2_sorted.bam"
    "SRR13890273_HISAT2_sorted.bam"
)

# Set the GTF annotation file
GTF_FILE="gencode.v37.annotation.gtf"

# Set the output directory
OUTPUT_DIR="htseq_counts_output"

# Make the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process each BAM file
for CURRENT_BAM in "${BAM_FILES[@]}"; do
    # Run HTSeq-Count
    htseq-count -f bam -r pos -s no -t exon -i gene_id "$CURRENT_BAM" "$GTF_FILE" > "${OUTPUT_DIR}/${CURRENT_BAM}_counts.txt"
done
