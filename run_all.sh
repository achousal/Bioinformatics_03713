#!/bin/bash
#SBATCH --job-name=run_all_submission
#SBATCH --output=logs/run_all_submission_%j.out
#SBATCH --error=logs/run_all_submission_%j.err
#SBATCH --time=96:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=9
#SBATCH --partition=RM-shared

# Create logs directory if missing
mkdir -p logs

echo "Submitting enhancer.sh..."
ENHANCER_JOBID=$(sbatch enhancer.sh | awk '{print $4}')
echo "Enhancer job submitted with Job ID: $ENHANCER_JOBID"

echo "Submitting map_halper.sh..."
MAP_HALPER_JOBID=$(sbatch map_halper.sh | awk '{print $4}')
echo "HALPER mapping job submitted with Job ID: $MAP_HALPER_JOBID"

# Wait for enhancer and mapping jobs to finish before submitting dependent jobs
echo "Waiting for enhancer and halper mapping jobs to complete before proceeding..."
DEPENDENCY="afterok:${ENHANCER_JOBID}:${MAP_HALPER_JOBID}"

# Submit cross_species.sh with dependency
echo "Submitting cross_species.sh with dependency on enhancer and mapping..."
# Replace paths with correct files!
CROSS_SPECIES_JOBID=$(sbatch --dependency=$DEPENDENCY --wrap="bash cross_species.sh \
  ./path/to/human_liver_peaks.bed \
  ./path/to/human_pancreas_peaks.bed \
  ./path/to/mouse_liver_peaks.bed \
  ./path/to/mouse_pancreas_peaks.bed \
  ./path/to/human_liver_to_mouse_halper.bed \
  ./path/to/human_pancreas_to_mouse_halper.bed \
  ./path/to/mouse_liver_to_human_halper.bed \
  ./path/to/mouse_pancreas_to_human_halper.bed" | awk '{print $4}')
echo "Cross-species job submitted with Job ID: $CROSS_SPECIES_JOBID"

# Submit meme-chip.sh with dependency
echo "Submitting meme-chip.sh with dependency on enhancer and mapping..."
MEME_CHIP_JOBID=$(sbatch --dependency=$DEPENDENCY meme-chip.sh | awk '{print $4}')
echo "MEME-ChIP job submitted with Job ID: $MEME_CHIP_JOBID"

echo "All jobs submitted!"
