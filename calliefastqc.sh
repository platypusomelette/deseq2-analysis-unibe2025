#!/bin/bash
#SBATCH --job-name=sample1blood
#SBATCH --output=/data/users/ksales/output_dir/fastqc_%j.out
#SBATCH --error=/data/users/ksales/output_dir/fastqc_%j.err
#SBATCH --time=02:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8

# Set paths
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"
INPUT_DIR="/data/users/ksales/reads_Blood"
OUTPUT_DIR="/data/users/ksales/output_dir"

# Run FastQC using Apptainer container
apptainer exec "$CONTAINER" fastqc \
    "$INPUT_DIR"/*.fastq* \
    -o "$OUTPUT_DIR" \
    -t $SLURM_CPUS_PER_TASK

echo "FastQC analysis complete!"
