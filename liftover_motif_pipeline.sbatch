#!/bin/bash
#SBATCH --job-name=liftover_motif_pipeline
#SBATCH --output=liftover_motif_%j.out
#SBATCH --error=liftover_motif_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=RM

# Load required modules and activate environment
module load anaconda3
conda activate bio

# --- User-specified input paths ---
read -p "Path to HAL alignment file: " HAL_FILE
read -p "Path to human peaks BED file: " HUMAN_PEAKS
read -p "Path to mouse genome FASTA (e.g., mm10.fa): " MOUSE_FASTA
read -p "Path to mouse open chromatin BED file: " MOUSE_OPEN_BED
read -p "Path to motif database (e.g., JASPAR .meme file): " MOTIF_DB
read -p "Output directory: " OUTPUT_DIR

# Derived output paths
LIFTOVER_BED="$OUTPUT_DIR/mouse_peaks_liftover.bed"
INTERSECT_BED="$OUTPUT_DIR/mouse_peaks_liftover_intersect_mouse_liver.bed"
FASTA_OUT="$OUTPUT_DIR/mouse_conserved_sequences.fa"
MEME_OUT="$OUTPUT_DIR/meme_output"
CONSERVATION_STATS="$OUTPUT_DIR/conservation_summary.txt"

mkdir -p "$OUTPUT_DIR"

# --- Step 1: Run halLiftover ---
echo "Running halLiftover..."
halLiftover "$HAL_FILE" \
  Human \
  "$HUMAN_PEAKS" \
  Mouse \
  "$LIFTOVER_BED"

# --- Step 2: Intersect with mouse liver peaks ---
echo "Running bedtools intersect..."
bedtools intersect -a "$LIFTOVER_BED" \
  -b "$MOUSE_OPEN_BED" \
  -u > "$INTERSECT_BED"

# --- Step 3: Extract FASTA sequences ---
echo "Extracting FASTA sequences..."
bedtools getfasta -fi "$MOUSE_FASTA" \
  -bed "$INTERSECT_BED" \
  -fo "$FASTA_OUT"


