#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=9
#SBATCH -J halper_mapping
#SBATCH -o halper_mapping_%j.out
#SBATCH -e halper_mapping_%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Set paths
HALPER_DIR="/ocean/projects/bio230007p/ssabata/repos/halLiftover-postprocessing"
OUTPUT_DIR="$PROJECT/mapped_peaks"
CACTUS_ALIGNMENT="/ocean/projects/bio230007p/ssabata/Alignments/10plusway-master.hal"

# Create output directory
mkdir -p $OUTPUT_DIR

# Activate conda environment with HALPER dependencies
module load anaconda3
source activate halper_env

# Map human liver peaks to mouse
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/ssabata/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
  -o $OUTPUT_DIR \
  -s human \
  -t mouse \
  -c $CACTUS_ALIGNMENT \
  -n human_liver

# Map human pancreas peaks to mouse
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/ssabata/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
  -o $OUTPUT_DIR \
  -s human \
  -t mouse \
  -c $CACTUS_ALIGNMENT \
  -n human_pancreas

# Map mouse liver peaks to human
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/ssabata/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
  -o $OUTPUT_DIR \
  -s mouse \
  -t human \
  -c $CACTUS_ALIGNMENT \
  -n mouse_liver

# Map mouse pancreas peaks to human
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/ssabata/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
  -o $OUTPUT_DIR \
  -s mouse \
  -t human \
  -c $CACTUS_ALIGNMENT \
  -n mouse_pancreas

echo "All mapping jobs completed"
