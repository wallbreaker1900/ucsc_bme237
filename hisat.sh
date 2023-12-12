#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --partition=256x44
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G 

input_files=("SRR13890278" "SRR13890277" "SRR13890275" "SRR13890274" "SRR13890272" "SRR13890271")
index_prefix="HISAT2_index"

# Loop through each input file and run hisat2
for input_file in "${input_files[@]}"; do
    input_fastq="${input_file}_trimmed.fastq"
    output_sam="${input_file}_HISAT2_output.sam"
    output_sorted_bam="${input_file}_HISAT2_sorted.bam"
    
    /hb/home/anbrooks/bme237_22SP/bin/HISAT2/hisat2-2.1.0/hisat2 -x "${index_prefix}" -U "${input_fastq}" -S "${output_sam}" -p 20
    
    samtools sort -@ 20 -o "${output_sorted_bam}" "${output_sam}"
    
    echo "Processed ${input_fastq}, saved output to ${output_sorted_bam}"
    
    rm "${output_sam}"
done