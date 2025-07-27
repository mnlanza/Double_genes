#!/bin/bash

# Exit on error, undefined var, or pipeline fail
set -euo pipefail

# === Usage Check ===
if [ $# -lt 3 ]; then
    echo "Usage: $0 JOB_ID INPUT_FASTA OUTPUT_DIR (full paths)"
    echo "Example: $0 double-genes double_gyrae/input/double_genes.fasta output_dir(full path)."
    exit 1
fi

job_id="$1"
input_fasta="$2"
output_dir="$3"

# Create required directories
mkdir -p ${output_dir}/output ${output_dir}/jobs ${output_dir}/figures

# Set version
job_version="v1"

# Create job directory before running anything
job_dir="${output_dir}/jobs/${job_id}-${job_version}"
mkdir -p "$job_dir"

# === Submit Job ===
echo "Submitting job..."
evo_gcp submit \
    --job "$job_id" \
    --output_type logits \
    --input_fasta "${input_fasta}" \
    --job_version "$job_version" \
    --wait

# === Download Results ===
echo "Downloading results..."
evo_gcp download \
    --job "$job_id" \
    --job_version "$job_version" \
    --jobs_dir "${output_dir}/jobs"

echo "Job completed successfully!"
echo "Results are available in: ${job_dir}"