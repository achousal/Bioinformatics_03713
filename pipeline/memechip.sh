#!/bin/bash
# memechip.sh - Run MEME-ChIP analysis on various peak subsets generated in Step 4.
#               Extracts sequences using bedtools getfasta and runs meme-chip.
#
# Usage: ./memechip.sh <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_genome_fa> <mouse_genome_fa> <motif_db_human> <motif_db_mouse> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>
#
# Inputs:
#   Arg 1:  <output_dir>                  - Base directory for pipeline results.
#   Arg 2:  <human_liver_peaks_gz>      - Gzipped narrowPeak: Original Human Liver peaks (used for 'full' set).
#   Arg 3:  <human_pancreas_peaks_gz>   - Gzipped narrowPeak: Original Human Pancreas peaks (used for 'full' set).
#   Arg 4:  <mouse_liver_peaks_gz>      - Gzipped narrowPeak: Original Mouse Liver peaks (used for 'full' set).
#   Arg 5:  <mouse_pancreas_peaks_gz>   - Gzipped narrowPeak: Original Mouse Pancreas peaks (used for 'full' set).
#   Arg 6:  <human_genome_fa>           - FASTA: Human reference genome.
#   Arg 7:  <mouse_genome_fa>           - FASTA: Mouse reference genome.
#   Arg 8:  <motif_db_human>            - MEME format: Human motif database.
#   Arg 9:  <motif_db_mouse>            - MEME format: Mouse motif database.
#   Arg 10: <mapped_hl_to_mm_bed>       - BED: Mapped Human Liver peaks (Mouse coords, from Step 1, for skip check).
#   Arg 11: <mapped_hp_to_mm_bed>       - BED: Mapped Human Pancreas peaks (Mouse coords, from Step 1, for skip check).
#   Arg 12: <mapped_ml_to_hg_bed>       - BED: Mapped Mouse Liver peaks (Human coords, from Step 1, for skip check).
#   Arg 13: <mapped_mp_to_hg_bed>       - BED: Mapped Mouse Pancreas peaks (Human coords, from Step 1, for skip check).
#   Implicit Dependencies (Expected files from Step 4 within <output_dir>/classification_and_comparisons/):
#     - classified_peaks/human/human_*.bed       (Promoters/Enhancers)
#     - classified_peaks/mouse/mouse_*.bed       (Promoters/Enhancers)
#     - region_comparisons/tissue_comparison/*   (Tissue shared/specific Enhancers)
#     - region_comparisons/species_comparison/*  (Species shared/specific Enhancers vs mapped)
#
# Outputs:
#   Directory: <output_dir>/classification_and_comparisons/full_peak_beds/ - Created by this script.
#   Files:     .../*.bed                                                   - BED3 versions of original peaks.
#   Directory: <output_dir>/sequences/                                     - Created by this script.
#   Files:     .../<subset_name>.fa                                        - FASTA sequences for each analyzed subset.
#   Directory: <output_dir>/meme_chip_results/                             - Created by this script.
#   Subdirs:   .../<subset_name>/                                          - MEME-ChIP output directories for each subset.
#   Files:     .../<subset_name>/meme-chip.html, meme_out/, etc.           - Standard MEME-ChIP results.
#   Files:     .../<subset_name>/meme-chip.log                             - Log file for each MEME-ChIP run.

set -e

# Parse command-line arguments
if [ $# -lt 13 ]; then
    echo "Usage: $0 <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_genome_fa> <mouse_genome_fa> <motif_db_human> <motif_db_mouse> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>"
    exit 1
fi

OUTPUT_DIR=$1
HUMAN_LIVER_PEAKS=$2 # Original narrowPeak.gz path (for full peak set)
HUMAN_PANCREAS_PEAKS=$3 # Original narrowPeak.gz path (for full peak set)
MOUSE_LIVER_PEAKS=$4 # Original narrowPeak.gz path (for full peak set)
MOUSE_PANCREAS_PEAKS=$5 # Original narrowPeak.gz path (for full peak set)
HUMAN_GENOME=$6
MOUSE_GENOME=$7
MOTIF_DB_HUMAN=$8
MOTIF_DB_MOUSE=$9
MAPPED_HL_TO_MM_PEAKS=${10} # Mapped BED path from Step 1 (for SKIP_SPECIES_STATS)
MAPPED_HP_TO_MM_PEAKS=${11} # Mapped BED path from Step 1 (for SKIP_SPECIES_STATS)
MAPPED_ML_TO_HG_PEAKS=${12} # Mapped BED path from Step 1 (for SKIP_SPECIES_STATS)
MAPPED_MP_TO_HG_PEAKS=${13} # Mapped BED path from Step 1 (for SKIP_SPECIES_STATS)

# Define directories based on main output directory
ANALYSIS_OUTPUT_DIR="$OUTPUT_DIR/classification_and_comparisons"
SEQUENCE_DIR="$OUTPUT_DIR/sequences"
MEME_RESULTS="$OUTPUT_DIR/meme_chip_results"
FULL_PEAKS_BED_DIR="$ANALYSIS_OUTPUT_DIR/full_peak_beds" # Dir for BED3 version of full peaks
CLASSIFIED_DIR="$ANALYSIS_OUTPUT_DIR/classified_peaks"
COMP_DIR="$ANALYSIS_OUTPUT_DIR/region_comparisons"

# Determine if species stats were skipped in Step 4 (based on existence of mapped files)
SKIP_SPECIES_STATS=0
if [[ ! -s "$MAPPED_ML_TO_HG_PEAKS" || ! -s "$MAPPED_MP_TO_HG_PEAKS" || ! -s "$MAPPED_HL_TO_MM_PEAKS" || ! -s "$MAPPED_HP_TO_MM_PEAKS" ]]; then
    echo "  Warning: One or more mapped peak BED files used for species comparison are missing or empty. Skipping MEME-ChIP on species-specific/shared sets derived from them." >&2
    SKIP_SPECIES_STATS=1
fi

# --- Reconstruct paths to Step 4 output BED files --- Needed for meme_inputs array
# Classified peaks
HUMAN_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_liver_promoters.bed"
HUMAN_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_liver_enhancers.bed"
HUMAN_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_pancreas_promoters.bed"
HUMAN_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_pancreas_enhancers.bed"
MOUSE_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_promoters.bed"
MOUSE_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_enhancers.bed"
MOUSE_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_promoters.bed"
MOUSE_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_enhancers.bed"
# Tissue comparison
HUMAN_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/human_enhancers_tissue_shared.bed"
HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_liver_enhancers_tissue_specific.bed"
HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_pancreas_enhancers_tissue_specific.bed"
MOUSE_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/mouse_enhancers_tissue_shared.bed"
MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_liver_enhancers_tissue_specific.bed"
MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_pancreas_enhancers_tissue_specific.bed"
# Species comparison (enhancers)
HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"
HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"


# --- The rest of the script remains largely the same, using the arguments and reconstructed paths --- 

# ===== STEP 5: Run MEME-ChIP Analysis =====
echo "===== STEP 5: Running MEME-ChIP motif analysis ====="

mkdir -p "$FULL_PEAKS_BED_DIR" "$SEQUENCE_DIR" "$MEME_RESULTS" # Ensure sequence/meme dirs exist

echo "  Preparing full peak set BED files (if needed)..."
# Create BED3 files for the full peak sets (category a) only if they don't exist
if [[ ! -s "$FULL_PEAKS_BED_DIR/human_liver_full.bed" ]]; then zcat "$HUMAN_LIVER_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/human_liver_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/human_pancreas_full.bed" ]]; then zcat "$HUMAN_PANCREAS_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/human_pancreas_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/mouse_liver_full.bed" ]]; then zcat "$MOUSE_LIVER_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/mouse_liver_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed" ]]; then zcat "$MOUSE_PANCREAS_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed"; fi

# Define the list of BED files for MEME-ChIP (categories a-g from Step 4 outputs)
declare -A meme_inputs=(
    # a. Full peak sets
    ["human_liver_full"]="$FULL_PEAKS_BED_DIR/human_liver_full.bed"
    ["human_pancreas_full"]="$FULL_PEAKS_BED_DIR/human_pancreas_full.bed"
    ["mouse_liver_full"]="$FULL_PEAKS_BED_DIR/mouse_liver_full.bed"
    ["mouse_pancreas_full"]="$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed"
    # b. Enhancers
    ["human_liver_enhancers"]="$HUMAN_LIVER_ENHANCERS_BED"
    ["human_pancreas_enhancers"]="$HUMAN_PANCREAS_ENHANCERS_BED"
    ["mouse_liver_enhancers"]="$MOUSE_LIVER_ENHANCERS_BED"
    ["mouse_pancreas_enhancers"]="$MOUSE_PANCREAS_ENHANCERS_BED"
    # c. Promoters
    ["human_liver_promoters"]="$HUMAN_LIVER_PROMOTERS_BED"
    ["human_pancreas_promoters"]="$HUMAN_PANCREAS_PROMOTERS_BED"
    ["mouse_liver_promoters"]="$MOUSE_LIVER_PROMOTERS_BED"
    ["mouse_pancreas_promoters"]="$MOUSE_PANCREAS_PROMOTERS_BED"
    # d. Enhancers shared across tissues
    ["human_enhancers_tissue_shared"]="$HUMAN_ENHANCERS_TISSUE_SHARED"
    ["mouse_enhancers_tissue_shared"]="$MOUSE_ENHANCERS_TISSUE_SHARED"
    # e. Enhancers specific to each tissue
    ["human_liver_enhancers_tissue_specific"]="$HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC"
    ["human_pancreas_enhancers_tissue_specific"]="$HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"
    ["mouse_liver_enhancers_tissue_specific"]="$MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC"
    ["mouse_pancreas_enhancers_tissue_specific"]="$MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"
    # f. Enhancers shared across species (vs mapped peaks)
    ["human_liver_enhancers_species_shared_hg"]="$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
    ["human_pancreas_enhancers_species_shared_hg"]="$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
    ["mouse_liver_enhancers_species_shared_mm"]="$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
    ["mouse_pancreas_enhancers_species_shared_mm"]="$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
    # g. Enhancers specific to each species (vs mapped peaks)
    ["human_liver_enhancers_species_specific_hg"]="$HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
    ["human_pancreas_enhancers_species_specific_hg"]="$HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
    ["mouse_liver_enhancers_species_specific_mm"]="$MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
    ["mouse_pancreas_enhancers_species_specific_mm"]="$MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
)

# Loop through the defined inputs and run MEME-ChIP
for name in "${!meme_inputs[@]}"; do
  BED_FILE="${meme_inputs[$name]}"
  FASTA_FILE="$SEQUENCE_DIR/${name}.fa"
  MEME_OUT_DIR="$MEME_RESULTS/${name}"

  # Determine species and genome/motif DB
  if [[ "$name" == human* ]]; then GENOME="$HUMAN_GENOME"; SPECIES="human"; MOTIF_DB="$MOTIF_DB_HUMAN";
  elif [[ "$name" == mouse* ]]; then GENOME="$MOUSE_GENOME"; SPECIES="mouse"; MOTIF_DB="$MOTIF_DB_MOUSE";
  else echo "    Warning: Cannot determine species for $name based on name. Skipping MEME-ChIP." >&2; continue; fi

  echo "    Processing $name for MEME-ChIP..."

  # Check if BED file exists and is not empty
  # For species comparison files, also check if the step was skipped
  if [[ "$SKIP_SPECIES_STATS" -eq 1 && "$name" == *species_* ]]; then echo "    Skipping $name: Species comparison step was skipped."; continue; fi
  if [[ ! -s "$BED_FILE" ]]; then echo "    Skipping $name: Input BED file $BED_FILE not found or empty."; continue; fi

  # Create FASTA file
  if [[ -s "$FASTA_FILE" ]]; then echo "      Skipping getfasta for $name: FASTA file exists."
  else echo "      Running getfasta for $name..."; bedtools getfasta -fi "$GENOME" -bed "$BED_FILE" -fo "$FASTA_FILE"; if [[ ! -s "$FASTA_FILE" ]]; then echo "    Skipping $name: Failed to create FASTA file or it is empty."; continue; fi; fi

  # Run MEME-ChIP
  if [[ -d "$MEME_OUT_DIR/meme_out" || -f "$MEME_OUT_DIR/meme-chip.html" ]]; then echo "      Skipping MEME-ChIP for $name: Output directory exists."
  else mkdir -p "$MEME_OUT_DIR"; echo "      Running MEME-ChIP for $name..."; meme-chip -oc "$MEME_OUT_DIR" -db "$MOTIF_DB" -meme-nmotifs 5 -meme-maxw 20 -seed 42 "$FASTA_FILE" &> "$MEME_OUT_DIR/meme-chip.log"; if [ $? -ne 0 ]; then echo "Error: MEME-ChIP failed for $name. Check log: $MEME_OUT_DIR/meme-chip.log" >&2; else echo "      Finished MEME-ChIP for $name."; fi; fi
done

echo "===== Completed Step 5 ====="