#(base) [villegmb@login509-02-r go_analysis]$ more creating_slurm_jobs_for_go_annotation.sh

#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J creating_slurm_jobs
#SBATCH -o creating_slurm_jobs.%J.out
#SBATCH -e creating_slurm_jobs.%J.err
#SBATCH --time=50:30:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24

# Define the directory containing blastp results
BLASTP_RESULT_DIR="/ibex/tmp/c2078/Heat_stress_analysis/scripts/go_analysis/filtered_blastp.outfmt6_trembl/tsv/merged_files/"

# Loop through the numbers 1 to 100
for i in {1..611}; do
    # Define the input blastp result file
    BLASTP_RESULT_FILE="${BLASTP_RESULT_DIR}merged_group_${i}.tsv"

    # Define the output folder
    OUTPUT_FOLDER="${BLASTP_RESULT_DIR}output_${i}_go_annotations/"

    # Create output folder if it doesn't exist
    mkdir -p "$OUTPUT_FOLDER"

    # Create SLURM job script
    cat > "running_annotation_${i}.sh" <<EOF
#!/bin/bash --login
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J running_annotation_${i}
#SBATCH -o running_annotation_${i}.%J.out
#SBATCH -e running_annotation_${i}.%J.err
#SBATCH --time=100:30:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24

# Define the input blastp result file
BLASTP_RESULT_FILE="${BLASTP_RESULT_FILE}"

# Define the input folder (same as the directory containing the blastp result file)
INPUT_FOLDER="$BLASTP_RESULT_DIR"

# Extract the name of the input folder
INPUT_FOLDER_NAME=\$(basename "\$BLASTP_RESULT_FILE" | cut -d '.' -f1)

# Define the output folder
OUTPUT_FOLDER="${OUTPUT_FOLDER}"

# Create output folder if it doesn't exist
mkdir -p "\$OUTPUT_FOLDER"

# Run the Python script
python create_go_annots_sprot_output.py ahem "\$BLASTP_RESULT_FILE" "\$OUTPUT_FOLDER" -n
EOF

    # Make the SLURM job script executable
    chmod +x "running_annotation_${i}.sh"
done
