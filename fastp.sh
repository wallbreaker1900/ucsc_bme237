#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=256x44 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G

input_files=("SRR13890278.fastq" "SRR13890277.fastq" "SRR13890275.fastq" "SRR13890274.fastq" "SRR13890272.fastq" "SRR13890271.fastq")

for input_file in "${input_files[@]}"; do
    output_file="${input_file%.fastq}_trimmed.fastq"
    
    fastp -i "$input_file" -o "$output_file" -w 20
    
    echo "Processed $input_file and saved to $output_file"
done
