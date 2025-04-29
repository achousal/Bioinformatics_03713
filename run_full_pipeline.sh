#!/bin/bash
#SBATCH --job-name=full_pipeline
#SBATCH --output=logs/full_pipeline_%j.out
#SBATCH --error=logs/full_pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=64GB

set -e

# Load necessary modules
module load bedtools
module load MEME-suite/5.4.1
module load samtools
module load anaconda3
conda init
source activate hal

export PATH=/jet/home/achousal/repos/hal/bin:${PATH}
export PYTHONPATH=/jet/home/achousal/repos/halLiftover-postprocessing:${PYTHONPATH}

# ===== CONFIGURATION =====
BASE_DIR="/ocean/projects/bio230007p/achousal"
HALPER_DIR="/jet/home/achousal/repos/halLiftover-postprocessing"
CACTUS_ALIGNMENT="/ocean/projects/bio230007p/achousal/ikaplow/Alignments/10plusway-master.hal"
CLASSIFIED_DIR="$BASE_DIR/promoters_enhancers"
GENOME_DIR="$BASE_DIR/ikaplow"
OUTPUT_DIR="$BASE_DIR/results"
MAPPED_DIR="$OUTPUT_DIR/mapped_peaks"
SEQUENCE_DIR="$OUTPUT_DIR/sequences"
MEME_RESULTS="$OUTPUT_DIR/meme_results"
mkdir -p logs "$OUTPUT_DIR" "$MAPPED_DIR" "$SEQUENCE_DIR" "$MEME_RESULTS"

# GTF files
HUMAN_GTF="$GENOME_DIR/HumanGenomeInfo/gencode.v47.annotation.gff3.gz"
MOUSE_GTF="$GENOME_DIR/MouseGenomeInfo/gencode.vM10.annotation.gff3.gz"

# Raw peak files
HUMAN_LIVER_PEAKS="$BASE_DIR/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
HUMAN_PANCREAS_PEAKS="$BASE_DIR/ikaplow/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_LIVER_PEAKS="$BASE_DIR/ikaplow/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_PANCREAS_PEAKS="$BASE_DIR/ikaplow/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"

# Genome FASTA files
HUMAN_GENOME="$GENOME_DIR/HumanGenomeInfo/hg38.fa"
MOUSE_GENOME="$GENOME_DIR/MouseGenomeInfo/mm10.fa"

# Motif databases
MOTIF_DB_HUMAN="$GENOME_DIR/CIS-BP_2.00/Homo_sapiens.meme"
MOTIF_DB_MOUSE="$GENOME_DIR/CIS-BP_2.00/Mus_musculus.meme"

# Index genomes
samtools faidx "$HUMAN_GENOME"
samtools faidx "$MOUSE_GENOME"

# ===== STEP 1: Map with HALPER =====
echo "Running HALPER liftover..."
mkdir -p "$MAPPED_DIR"

function run_halper() {
  local species_source=$1
  local species_target=$2
  local input_peaks=$3
  local name=$4

  local tmp_bed="/tmp/${name}_peaks.bed"
  zcat "$input_peaks" > "$tmp_bed"

  bash "$HALPER_DIR/halper_map_peak_orthologs.sh" \
    -b "$tmp_bed" \
    -o "$MAPPED_DIR" \
    -s "$species_source" \
    -t "$species_target" \
    -c "$CACTUS_ALIGNMENT" \
    -n "$name"

  zcat "$MAPPED_DIR/${name}.${species_source}To${species_target}.HALPER.narrowPeak.gz" | cut -f1-3 > "$MAPPED_DIR/${name}_${species_target}_mapped.bed"
}

run_halper Human Mouse "$HUMAN_LIVER_PEAKS" human_liver
run_halper Human Mouse "$HUMAN_PANCREAS_PEAKS" human_pancreas
run_halper Mouse Human "$MOUSE_LIVER_PEAKS" mouse_liver
run_halper Mouse Human "$MOUSE_PANCREAS_PEAKS" mouse_pancreas

# ===== CAPTURE HALPER OUTPUT PATHS =====
HUMAN_LIVER_TO_MOUSE_HALPER="$MAPPED_DIR/human_liver.HumanToMouse.HALPER.narrowPeak.gz"
HUMAN_PANCREAS_TO_MOUSE_HALPER="$MAPPED_DIR/human_pancreas.HumanToMouse.HALPER.narrowPeak.gz"
MOUSE_LIVER_TO_HUMAN_HALPER="$MAPPED_DIR/mouse_liver.MouseToHuman.HALPER.narrowPeak.gz"
MOUSE_PANCREAS_TO_HUMAN_HALPER="$MAPPED_DIR/mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz"

# ===== STEP 2: Cross-species and Tissue Comparison =====
echo "Running tissue-specific and conservation analyses..."

# ===== CONFIGURATION =====
# Load required modules
module load bedtools/2.30.0 

# Create output directories
RESULTS_DIR="./results"
CROSS_SPECIES_DIR="$RESULTS_DIR/cross_species"
TISSUE_COMP_DIR="$RESULTS_DIR/tissue_comparison"
PROCESSED_DIR="$CROSS_SPECIES_DIR/processed_files"
mkdir -p $RESULTS_DIR $CROSS_SPECIES_DIR $TISSUE_COMP_DIR $PROCESSED_DIR


# Log input files
echo "===== INPUT FILES ====="
echo "Human liver peaks: $HUMAN_LIVER_PEAKS"
echo "Human pancreas peaks: $HUMAN_PANCREAS_PEAKS"
echo "Mouse liver peaks: $MOUSE_LIVER_PEAKS"
echo "Mouse pancreas peaks: $MOUSE_PANCREAS_PEAKS"
echo "Human liver to mouse HALPER: $HUMAN_LIVER_TO_MOUSE_HALPER"
echo "Human pancreas to mouse HALPER: $HUMAN_PANCREAS_TO_MOUSE_HALPER"
echo "Mouse liver to human HALPER: $MOUSE_LIVER_TO_HUMAN_HALPER"
echo "Mouse pancreas to human HALPER: $MOUSE_PANCREAS_TO_HUMAN_HALPER"
echo ""

# Check if input files exist
for file in "$HUMAN_LIVER_PEAKS" "$HUMAN_PANCREAS_PEAKS" "$MOUSE_LIVER_PEAKS" "$MOUSE_PANCREAS_PEAKS" "$HUMAN_LIVER_TO_MOUSE_HALPER" "$HUMAN_PANCREAS_TO_MOUSE_HALPER" "$MOUSE_LIVER_TO_HUMAN_HALPER" "$MOUSE_PANCREAS_TO_HUMAN_HALPER"; do
    if [ ! -f "$file" ]; then
        echo "Error: Input file $file does not exist"
        exit 1
    fi
done

# Function to safely calculate percentages
calc_pct() {
    local num=$1
    local denom=$2
    if [ "$denom" -eq 0 ]; then
        echo "0.00"
    else
        echo "scale=2; $num / $denom * 100" | bc
    fi
}

# ===== STEP 1: WITHIN-SPECIES TISSUE COMPARISON =====
echo "===== STARTING TISSUE COMPARISON ANALYSIS ====="

# Convert narrowPeak files to BED format
echo "Converting narrowPeak files to BED format..."
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$HUMAN_LIVER_PEAKS" > "$TISSUE_COMP_DIR/human_liver.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$HUMAN_PANCREAS_PEAKS" > "$TISSUE_COMP_DIR/human_pancreas.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$MOUSE_LIVER_PEAKS" > "$TISSUE_COMP_DIR/mouse_liver.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$MOUSE_PANCREAS_PEAKS" > "$TISSUE_COMP_DIR/mouse_pancreas.bed"

# Human tissue comparison
echo "Finding shared and tissue-specific regions in human..."
HUMAN_SHARED="$TISSUE_COMP_DIR/human_shared_between_liver_and_pancreas.bed"
HUMAN_LIVER_SPECIFIC="$TISSUE_COMP_DIR/human_liver_specific.bed"
HUMAN_PANCREAS_SPECIFIC="$TISSUE_COMP_DIR/human_pancreas_specific.bed"

bedtools intersect -a "$TISSUE_COMP_DIR/human_liver.bed" -b "$TISSUE_COMP_DIR/human_pancreas.bed" -u > "$HUMAN_SHARED"
bedtools intersect -a "$TISSUE_COMP_DIR/human_liver.bed" -b "$TISSUE_COMP_DIR/human_pancreas.bed" -v > "$HUMAN_LIVER_SPECIFIC"
bedtools intersect -a "$TISSUE_COMP_DIR/human_pancreas.bed" -b "$TISSUE_COMP_DIR/human_liver.bed" -v > "$HUMAN_PANCREAS_SPECIFIC"

# Mouse tissue comparison
echo "Finding shared and tissue-specific regions in mouse..."
MOUSE_SHARED="$TISSUE_COMP_DIR/mouse_shared_between_liver_and_pancreas.bed"
MOUSE_LIVER_SPECIFIC="$TISSUE_COMP_DIR/mouse_liver_specific.bed"
MOUSE_PANCREAS_SPECIFIC="$TISSUE_COMP_DIR/mouse_pancreas_specific.bed"

bedtools intersect -a "$TISSUE_COMP_DIR/mouse_liver.bed" -b "$TISSUE_COMP_DIR/mouse_pancreas.bed" -u > "$MOUSE_SHARED"
bedtools intersect -a "$TISSUE_COMP_DIR/mouse_liver.bed" -b "$TISSUE_COMP_DIR/mouse_pancreas.bed" -v > "$MOUSE_LIVER_SPECIFIC"
bedtools intersect -a "$TISSUE_COMP_DIR/mouse_pancreas.bed" -b "$TISSUE_COMP_DIR/mouse_liver.bed" -v > "$MOUSE_PANCREAS_SPECIFIC"

# Calculate tissue sharing statistics
echo "Calculating tissue sharing statistics..."

# Count regions
HUMAN_LIVER_COUNT=$(wc -l < "$TISSUE_COMP_DIR/human_liver.bed")
HUMAN_PANCREAS_COUNT=$(wc -l < "$TISSUE_COMP_DIR/human_pancreas.bed")
HUMAN_SHARED_COUNT=$(wc -l < "$HUMAN_SHARED")
HUMAN_LIVER_SPECIFIC_COUNT=$(wc -l < "$HUMAN_LIVER_SPECIFIC")
HUMAN_PANCREAS_SPECIFIC_COUNT=$(wc -l < "$HUMAN_PANCREAS_SPECIFIC")

MOUSE_LIVER_COUNT=$(wc -l < "$TISSUE_COMP_DIR/mouse_liver.bed")
MOUSE_PANCREAS_COUNT=$(wc -l < "$TISSUE_COMP_DIR/mouse_pancreas.bed")
MOUSE_SHARED_COUNT=$(wc -l < "$MOUSE_SHARED")
MOUSE_LIVER_SPECIFIC_COUNT=$(wc -l < "$MOUSE_LIVER_SPECIFIC")
MOUSE_PANCREAS_SPECIFIC_COUNT=$(wc -l < "$MOUSE_PANCREAS_SPECIFIC")

# Calculate percentages
HUMAN_LIVER_SHARED_PCT=$(calc_pct $HUMAN_SHARED_COUNT $HUMAN_LIVER_COUNT)
HUMAN_PANCREAS_SHARED_PCT=$(calc_pct $HUMAN_SHARED_COUNT $HUMAN_PANCREAS_COUNT)
HUMAN_LIVER_SPECIFIC_PCT=$(calc_pct $HUMAN_LIVER_SPECIFIC_COUNT $HUMAN_LIVER_COUNT)
HUMAN_PANCREAS_SPECIFIC_PCT=$(calc_pct $HUMAN_PANCREAS_SPECIFIC_COUNT $HUMAN_PANCREAS_COUNT)

MOUSE_LIVER_SHARED_PCT=$(calc_pct $MOUSE_SHARED_COUNT $MOUSE_LIVER_COUNT)
MOUSE_PANCREAS_SHARED_PCT=$(calc_pct $MOUSE_SHARED_COUNT $MOUSE_PANCREAS_COUNT)
MOUSE_LIVER_SPECIFIC_PCT=$(calc_pct $MOUSE_LIVER_SPECIFIC_COUNT $MOUSE_LIVER_COUNT)
MOUSE_PANCREAS_SPECIFIC_PCT=$(calc_pct $MOUSE_PANCREAS_SPECIFIC_COUNT $MOUSE_PANCREAS_COUNT)

# Create tissue comparison summary - AVOIDING PARENTHESES
TISSUE_SUMMARY_FILE="$TISSUE_COMP_DIR/tissue_comparison_summary.txt"

echo "===== TISSUE COMPARISON SUMMARY =====" > "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "HUMAN TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"
echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $HUMAN_LIVER_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Total pancreas regions: $HUMAN_PANCREAS_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $HUMAN_SHARED_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $HUMAN_LIVER_SPECIFIC_COUNT - ${HUMAN_LIVER_SPECIFIC_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas-specific regions: $HUMAN_PANCREAS_SPECIFIC_COUNT - ${HUMAN_PANCREAS_SPECIFIC_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: ${HUMAN_LIVER_SHARED_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas sharing with liver: ${HUMAN_PANCREAS_SHARED_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "MOUSE TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"
echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $MOUSE_LIVER_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Total pancreas regions: $MOUSE_PANCREAS_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $MOUSE_SHARED_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $MOUSE_LIVER_SPECIFIC_COUNT - ${MOUSE_LIVER_SPECIFIC_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas-specific regions: $MOUSE_PANCREAS_SPECIFIC_COUNT - ${MOUSE_PANCREAS_SPECIFIC_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: ${MOUSE_LIVER_SHARED_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas sharing with liver: ${MOUSE_PANCREAS_SHARED_PCT}%" >> "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "Tissue comparison complete. Results in $TISSUE_COMP_DIR"

# ===== PART 1 COMPLETE =====
# ===== STEP 2: CROSS-SPECIES CONSERVATION ANALYSIS =====
echo "===== STARTING CROSS-SPECIES CONSERVATION ANALYSIS ====="

# Extract properly formatted coordinates from HALPER output
echo "Processing HALPER output files..."

# Function to process HALPER files and extract coordinates
process_halper_file() {
    local input_file=$1
    local output_file=$2
    
    if [[ "$input_file" == *.gz ]]; then
        # If file is gzipped
        zcat "$input_file" | awk '{
            split($4, coords, ":");
            if (length(coords) >= 3) {
                split(coords[2], range, "-");
                if (length(range) >= 2) {
                    print coords[1]"\t"range[1]"\t"range[2]"\t"$1":"$2"-"$3":"coords[3]"\t"$10
                }
            }
        }' > "$output_file"
    else
        # If file is not gzipped
        awk '{
            split($4, coords, ":");
            if (length(coords) >= 3) {
                split(coords[2], range, "-");
                if (length(range) >= 2) {
                    print coords[1]"\t"range[1]"\t"range[2]"\t"$1":"$2"-"$3":"coords[3]"\t"$10
                }
            }
        }' "$input_file" > "$output_file"
    fi
    
    # Count lines and print info
    line_count=$(wc -l < "$output_file")
    echo "  Processed $input_file: $line_count regions"
}

# Process HALPER files
process_halper_file "$HUMAN_LIVER_TO_MOUSE_HALPER" "$PROCESSED_DIR/human_liver_to_mouse.bed"
process_halper_file "$HUMAN_PANCREAS_TO_MOUSE_HALPER" "$PROCESSED_DIR/human_pancreas_to_mouse.bed"
process_halper_file "$MOUSE_LIVER_TO_HUMAN_HALPER" "$PROCESSED_DIR/mouse_liver_to_human.bed"
process_halper_file "$MOUSE_PANCREAS_TO_HUMAN_HALPER" "$PROCESSED_DIR/mouse_pancreas_to_human.bed"

# Find conserved regions using the processed HALPER output
echo "Finding conserved regions across species..."

# Human liver regions conserved in mouse
HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_liver.bed"
bedtools intersect -a "$PROCESSED_DIR/human_liver_to_mouse.bed" -b "$TISSUE_COMP_DIR/mouse_liver.bed" -u > "$HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER"

HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_pancreas.bed"
bedtools intersect -a "$PROCESSED_DIR/human_liver_to_mouse.bed" -b "$TISSUE_COMP_DIR/mouse_pancreas.bed" -u > "$HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS"

# Human pancreas regions conserved in mouse
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_pancreas.bed"
bedtools intersect -a "$PROCESSED_DIR/human_pancreas_to_mouse.bed" -b "$TISSUE_COMP_DIR/mouse_pancreas.bed" -u > "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS"

HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_liver.bed"
bedtools intersect -a "$PROCESSED_DIR/human_pancreas_to_mouse.bed" -b "$TISSUE_COMP_DIR/mouse_liver.bed" -u > "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER"

# Mouse liver regions conserved in human
MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER="$CROSS_SPECIES_DIR/mouse_liver_conserved_in_human_liver.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_liver_to_human.bed" -b "$TISSUE_COMP_DIR/human_liver.bed" -u > "$MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER"

MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS="$CROSS_SPECIES_DIR/mouse_liver_conserved_in_human_pancreas.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_liver_to_human.bed" -b "$TISSUE_COMP_DIR/human_pancreas.bed" -u > "$MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS"

# Mouse pancreas regions conserved in human
MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS="$CROSS_SPECIES_DIR/mouse_pancreas_conserved_in_human_pancreas.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_pancreas_to_human.bed" -b "$TISSUE_COMP_DIR/human_pancreas.bed" -u > "$MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS"

MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER="$CROSS_SPECIES_DIR/mouse_pancreas_conserved_in_human_liver.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_pancreas_to_human.bed" -b "$TISSUE_COMP_DIR/human_liver.bed" -u > "$MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER"

# Find non-conserved/unique regions
echo "Finding species-specific unique regions..."

# Create combined bed files for all tissues in each species
cat "$TISSUE_COMP_DIR/mouse_liver.bed" "$TISSUE_COMP_DIR/mouse_pancreas.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$PROCESSED_DIR/mouse_all_tissues.bed"
cat "$TISSUE_COMP_DIR/human_liver.bed" "$TISSUE_COMP_DIR/human_pancreas.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$PROCESSED_DIR/human_all_tissues.bed"

# Human regions not conserved in mouse
HUMAN_LIVER_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_liver_not_conserved_in_mouse.bed"
bedtools intersect -a "$PROCESSED_DIR/human_liver_to_mouse.bed" -b "$PROCESSED_DIR/mouse_all_tissues.bed" -v > "$HUMAN_LIVER_NOT_CONSERVED"

HUMAN_PANCREAS_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_pancreas_not_conserved_in_mouse.bed"
bedtools intersect -a "$PROCESSED_DIR/human_pancreas_to_mouse.bed" -b "$PROCESSED_DIR/mouse_all_tissues.bed" -v > "$HUMAN_PANCREAS_NOT_CONSERVED"

# Mouse regions not conserved in human
MOUSE_LIVER_NOT_CONSERVED="$CROSS_SPECIES_DIR/mouse_liver_not_conserved_in_human.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_liver_to_human.bed" -b "$PROCESSED_DIR/human_all_tissues.bed" -v > "$MOUSE_LIVER_NOT_CONSERVED"

MOUSE_PANCREAS_NOT_CONSERVED="$CROSS_SPECIES_DIR/mouse_pancreas_not_conserved_in_human.bed"
bedtools intersect -a "$PROCESSED_DIR/mouse_pancreas_to_human.bed" -b "$PROCESSED_DIR/human_all_tissues.bed" -v > "$MOUSE_PANCREAS_NOT_CONSERVED"

# Calculate cross-species conservation statistics
echo "Calculating cross-species conservation statistics..."

# Count regions for cross-species analysis
HUMAN_LIVER_TO_MOUSE_COUNT=$(wc -l < "$PROCESSED_DIR/human_liver_to_mouse.bed")
HUMAN_PANCREAS_TO_MOUSE_COUNT=$(wc -l < "$PROCESSED_DIR/human_pancreas_to_mouse.bed")
MOUSE_LIVER_TO_HUMAN_COUNT=$(wc -l < "$PROCESSED_DIR/mouse_liver_to_human.bed")
MOUSE_PANCREAS_TO_HUMAN_COUNT=$(wc -l < "$PROCESSED_DIR/mouse_pancreas_to_human.bed")

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER")
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER")

MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER_COUNT=$(wc -l < "$MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER")
MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS_COUNT=$(wc -l < "$MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS")
MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS_COUNT=$(wc -l < "$MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS")
MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER_COUNT=$(wc -l < "$MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER")

HUMAN_LIVER_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_LIVER_NOT_CONSERVED")
HUMAN_PANCREAS_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_PANCREAS_NOT_CONSERVED")
MOUSE_LIVER_NOT_CONSERVED_COUNT=$(wc -l < "$MOUSE_LIVER_NOT_CONSERVED")
MOUSE_PANCREAS_NOT_CONSERVED_COUNT=$(wc -l < "$MOUSE_PANCREAS_NOT_CONSERVED")

# Calculate percentages
HUMAN_LIVER_LIFTOVER_RATE=$(calc_pct $HUMAN_LIVER_TO_MOUSE_COUNT $HUMAN_LIVER_COUNT)
HUMAN_PANCREAS_LIFTOVER_RATE=$(calc_pct $HUMAN_PANCREAS_TO_MOUSE_COUNT $HUMAN_PANCREAS_COUNT)
MOUSE_LIVER_LIFTOVER_RATE=$(calc_pct $MOUSE_LIVER_TO_HUMAN_COUNT $MOUSE_LIVER_COUNT)
MOUSE_PANCREAS_LIFTOVER_RATE=$(calc_pct $MOUSE_PANCREAS_TO_HUMAN_COUNT $MOUSE_PANCREAS_COUNT)

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE=$(calc_pct $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT $HUMAN_LIVER_TO_MOUSE_COUNT)
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(calc_pct $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT $HUMAN_LIVER_TO_MOUSE_COUNT)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(calc_pct $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT $HUMAN_PANCREAS_TO_MOUSE_COUNT)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE=$(calc_pct $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT $HUMAN_PANCREAS_TO_MOUSE_COUNT)

MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER_RATE=$(calc_pct $MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER_COUNT $MOUSE_LIVER_TO_HUMAN_COUNT)
MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS_RATE=$(calc_pct $MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS_COUNT $MOUSE_LIVER_TO_HUMAN_COUNT)
MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS_RATE=$(calc_pct $MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS_COUNT $MOUSE_PANCREAS_TO_HUMAN_COUNT)
MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER_RATE=$(calc_pct $MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER_COUNT $MOUSE_PANCREAS_TO_HUMAN_COUNT)

HUMAN_LIVER_UNIQUE_RATE=$(calc_pct $HUMAN_LIVER_NOT_CONSERVED_COUNT $HUMAN_LIVER_TO_MOUSE_COUNT)
HUMAN_PANCREAS_UNIQUE_RATE=$(calc_pct $HUMAN_PANCREAS_NOT_CONSERVED_COUNT $HUMAN_PANCREAS_TO_MOUSE_COUNT)
MOUSE_LIVER_UNIQUE_RATE=$(calc_pct $MOUSE_LIVER_NOT_CONSERVED_COUNT $MOUSE_LIVER_TO_HUMAN_COUNT)
MOUSE_PANCREAS_UNIQUE_RATE=$(calc_pct $MOUSE_PANCREAS_NOT_CONSERVED_COUNT $MOUSE_PANCREAS_TO_HUMAN_COUNT)

# Create cross-species conservation summary - AVOIDING PARENTHESES
CROSS_SPECIES_SUMMARY_FILE="$CROSS_SPECIES_DIR/cross_species_conservation_summary.txt"

echo "===== CROSS-SPECIES CONSERVATION SUMMARY =====" > "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "HUMAN TO MOUSE CONSERVATION STATISTICS" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total human liver regions: $HUMAN_LIVER_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total human pancreas regions: $HUMAN_PANCREAS_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions mapped to mouse: $HUMAN_LIVER_TO_MOUSE_COUNT - ${HUMAN_LIVER_LIFTOVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions mapped to mouse: $HUMAN_PANCREAS_TO_MOUSE_COUNT - ${HUMAN_PANCREAS_LIFTOVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "HUMAN TO MOUSE CONSERVATION RATES" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in mouse liver: $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT - ${HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in mouse pancreas: $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT - ${HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse pancreas: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT - ${HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse liver: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT - ${HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "MOUSE TO HUMAN CONSERVATION STATISTICS" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total mouse liver regions: $MOUSE_LIVER_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total mouse pancreas regions: $MOUSE_PANCREAS_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver regions mapped to human: $MOUSE_LIVER_TO_HUMAN_COUNT - ${MOUSE_LIVER_LIFTOVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse pancreas regions mapped to human: $MOUSE_PANCREAS_TO_HUMAN_COUNT - ${MOUSE_PANCREAS_LIFTOVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "MOUSE TO HUMAN CONSERVATION RATES" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver regions conserved in human liver: $MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER_COUNT - ${MOUSE_LIVER_CONSERVED_IN_HUMAN_LIVER_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver regions conserved in human pancreas: $MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS_COUNT - ${MOUSE_LIVER_CONSERVED_IN_HUMAN_PANCREAS_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse pancreas regions conserved in human pancreas: $MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS_COUNT - ${MOUSE_PANCREAS_CONSERVED_IN_HUMAN_PANCREAS_RATE}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse pancreas regions conserved in human liver: $MOUSE_PANCREAS_CONSERVED_IN_HUMAN_LIVER_COUNT - ${

# ===== STEP 3: Classify Enhancers vs Promoters =====
echo "Classifying peaks as promoters or enhancers..."
mkdir -p "$CLASSIFIED_DIR/human" "$CLASSIFIED_DIR/mouse"

# Generate promoter BEDs
zcat "$HUMAN_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$CLASSIFIED_DIR/human_promoters.bed"
zcat "$MOUSE_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$CLASSIFIED_DIR/mouse_promoters.bed"

# Classify peaks
function classify_peaks() {
  local input_peak=$1
  local promoter_bed=$2
  local output_dir=$3
  local name=$4

  bedtools intersect -a "$input_peak" -b "$promoter_bed" -u > "$output_dir/${name}_promoters.bed"
  bedtools intersect -a "$input_peak" -b "$promoter_bed" -v > "$output_dir/${name}_enhancers.bed"
}

classify_peaks "$HUMAN_LIVER_PEAKS" "$CLASSIFIED_DIR/human_promoters.bed" "$CLASSIFIED_DIR/human" human_liver
classify_peaks "$HUMAN_PANCREAS_PEAKS" "$CLASSIFIED_DIR/human_promoters.bed" "$CLASSIFIED_DIR/human" human_pancreas
classify_peaks "$MOUSE_LIVER_PEAKS" "$CLASSIFIED_DIR/mouse_promoters.bed" "$CLASSIFIED_DIR/mouse" mouse_liver
classify_peaks "$MOUSE_PANCREAS_PEAKS" "$CLASSIFIED_DIR/mouse_promoters.bed" "$CLASSIFIED_DIR/mouse" mouse_pancreas

# ===== STEP 4: Extract Sequences and Run MEME-ChIP =====
echo "Running MEME-ChIP motif analysis..."

samples=(
  "human liver promoters"
  "human liver enhancers"
  "human pancreas promoters"
  "human pancreas enhancers"
  "mouse liver promoters"
  "mouse liver enhancers"
  "mouse pancreas promoters"
  "mouse pancreas enhancers"
)

for sample in "${samples[@]}"; do
  read -r SPECIES TISSUE REGION <<< "$sample"

  BED_FILE="$CLASSIFIED_DIR/${SPECIES}/${SPECIES}_${TISSUE}_${REGION}.bed"
  FASTA_FILE="$SEQUENCE_DIR/${SPECIES}_${TISSUE}_${REGION}.fa"
  OUTPUT_DIR="$MEME_RESULTS/${SPECIES}_${TISSUE}_${REGION}"

  mkdir -p "$OUTPUT_DIR"

  if [[ "$SPECIES" == "human" ]]; then
    GENOME="$HUMAN_GENOME"
    MOTIF_DB="$MOTIF_DB_HUMAN"
  else
    GENOME="$MOUSE_GENOME"
    MOTIF_DB="$MOTIF_DB_MOUSE"
  fi

  bedtools getfasta -fi "$GENOME" -bed "$BED_FILE" -fo "$FASTA_FILE"

  if [[ ! -s "$FASTA_FILE" ]]; then
    echo "Skipping $SPECIES $TISSUE $REGION: empty FASTA."
    continue
  fi

  meme-chip \
    -oc "$OUTPUT_DIR" \
    -db "$MOTIF_DB" \
    -meme-nmotifs 5 \
    -maxw 20 \
    "$FASTA_FILE"

  echo "Finished motif discovery for $sample."
done

# ===== DONE =====
echo "Pipeline completed successfully."
