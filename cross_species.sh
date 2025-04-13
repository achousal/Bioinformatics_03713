#!/bin/bash
#SBATCH --job-name=cross_species_analysis
#SBATCH --output=cross_species_%j.out
#SBATCH --error=cross_species_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=RM

# ===== CONFIGURATION =====
# Load required modules
module load bedtools/2.30.0 

# Create output directories
RESULTS_DIR="./results"
CROSS_SPECIES_DIR="$RESULTS_DIR/cross_species"
TISSUE_COMP_DIR="$RESULTS_DIR/tissue_comparison"
PROCESSED_DIR="$CROSS_SPECIES_DIR/processed_files"
mkdir -p $RESULTS_DIR $CROSS_SPECIES_DIR $TISSUE_COMP_DIR $PROCESSED_DIR

# Define input files
# Original peak files (change these paths to match your file locations)
HUMAN_LIVER_PEAKS="narrowPeak/human_liver.narrowPeak"
HUMAN_PANCREAS_PEAKS="narrowPeak/human_pancreas.narrowPeak"
MOUSE_LIVER_PEAKS="narrowPeak/mouse_liver.narrowPeak"
MOUSE_PANCREAS_PEAKS="narrowPeak/mouse_pancreas.narrowPeak"

# HALPER output files (change these paths to match your file locations)
HUMAN_LIVER_TO_MOUSE_HALPER="halper_result/idr.optimal_peak_liver.HumanToMouse.HALPER.narrowPeak.gz"
HUMAN_PANCREAS_TO_MOUSE_HALPER="halper_result/idr.optimal_peak_pancreas.HumanToMouse.HALPER.narrowPeak.gz"

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
HUMAN_LIVER_SHARED_PCT=$(echo "scale=2; $HUMAN_SHARED_COUNT / $HUMAN_LIVER_COUNT * 100" | bc)
HUMAN_PANCREAS_SHARED_PCT=$(echo "scale=2; $HUMAN_SHARED_COUNT / $HUMAN_PANCREAS_COUNT * 100" | bc)
HUMAN_LIVER_SPECIFIC_PCT=$(echo "scale=2; $HUMAN_LIVER_SPECIFIC_COUNT / $HUMAN_LIVER_COUNT * 100" | bc)
HUMAN_PANCREAS_SPECIFIC_PCT=$(echo "scale=2; $HUMAN_PANCREAS_SPECIFIC_COUNT / $HUMAN_PANCREAS_COUNT * 100" | bc)

MOUSE_LIVER_SHARED_PCT=$(echo "scale=2; $MOUSE_SHARED_COUNT / $MOUSE_LIVER_COUNT * 100" | bc)
MOUSE_PANCREAS_SHARED_PCT=$(echo "scale=2; $MOUSE_SHARED_COUNT / $MOUSE_PANCREAS_COUNT * 100" | bc)
MOUSE_LIVER_SPECIFIC_PCT=$(echo "scale=2; $MOUSE_LIVER_SPECIFIC_COUNT / $MOUSE_LIVER_COUNT * 100" | bc)
MOUSE_PANCREAS_SPECIFIC_PCT=$(echo "scale=2; $MOUSE_PANCREAS_SPECIFIC_COUNT / $MOUSE_PANCREAS_COUNT * 100" | bc)

# Create tissue comparison summary
TISSUE_SUMMARY_FILE="$TISSUE_COMP_DIR/tissue_comparison_summary.txt"

echo "===== TISSUE COMPARISON SUMMARY =====" > "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "HUMAN TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"
echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $HUMAN_LIVER_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Total pancreas regions: $HUMAN_PANCREAS_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $HUMAN_SHARED_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $HUMAN_LIVER_SPECIFIC_COUNT ($HUMAN_LIVER_SPECIFIC_PCT%)" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas-specific regions: $HUMAN_PANCREAS_SPECIFIC_COUNT ($HUMAN_PANCREAS_SPECIFIC_PCT%)" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: $HUMAN_LIVER_SHARED_PCT%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas sharing with liver: $HUMAN_PANCREAS_SHARED_PCT%" >> "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "MOUSE TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"
echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $MOUSE_LIVER_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Total pancreas regions: $MOUSE_PANCREAS_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $MOUSE_SHARED_COUNT" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $MOUSE_LIVER_SPECIFIC_COUNT ($MOUSE_LIVER_SPECIFIC_PCT%)" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas-specific regions: $MOUSE_PANCREAS_SPECIFIC_COUNT ($MOUSE_PANCREAS_SPECIFIC_PCT%)" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: $MOUSE_LIVER_SHARED_PCT%" >> "$TISSUE_SUMMARY_FILE"
echo "Pancreas sharing with liver: $MOUSE_PANCREAS_SHARED_PCT%" >> "$TISSUE_SUMMARY_FILE"
echo "" >> "$TISSUE_SUMMARY_FILE"

echo "Tissue comparison complete. Results in $TISSUE_COMP_DIR"

# ===== STEP 2: CROSS-SPECIES CONSERVATION ANALYSIS =====
echo "===== STARTING CROSS-SPECIES CONSERVATION ANALYSIS ====="

# Extract properly formatted mouse coordinates from HALPER output
echo "Processing HALPER output files..."

# Function to process HALPER files and extract mouse coordinates
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

# Find non-conserved regions
HUMAN_LIVER_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_liver_not_conserved_in_mouse.bed"
cat "$TISSUE_COMP_DIR/mouse_liver.bed" "$TISSUE_COMP_DIR/mouse_pancreas.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$PROCESSED_DIR/mouse_all_tissues.bed"
bedtools intersect -a "$PROCESSED_DIR/human_liver_to_mouse.bed" -b "$PROCESSED_DIR/mouse_all_tissues.bed" -v > "$HUMAN_LIVER_NOT_CONSERVED"

HUMAN_PANCREAS_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_pancreas_not_conserved_in_mouse.bed"
bedtools intersect -a "$PROCESSED_DIR/human_pancreas_to_mouse.bed" -b "$PROCESSED_DIR/mouse_all_tissues.bed" -v > "$HUMAN_PANCREAS_NOT_CONSERVED"

# Calculate cross-species conservation statistics
echo "Calculating cross-species conservation statistics..."

# Count regions for cross-species analysis
HUMAN_LIVER_COUNT=$(wc -l < "$TISSUE_COMP_DIR/human_liver.bed")
HUMAN_PANCREAS_COUNT=$(wc -l < "$TISSUE_COMP_DIR/human_pancreas.bed")
HUMAN_LIVER_TO_MOUSE_COUNT=$(wc -l < "$PROCESSED_DIR/human_liver_to_mouse.bed")
HUMAN_PANCREAS_TO_MOUSE_COUNT=$(wc -l < "$PROCESSED_DIR/human_pancreas_to_mouse.bed")

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER")
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER")

HUMAN_LIVER_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_LIVER_NOT_CONSERVED")
HUMAN_PANCREAS_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_PANCREAS_NOT_CONSERVED")

# Calculate percentages
HUMAN_LIVER_LIFTOVER_RATE=$(echo "scale=2; $HUMAN_LIVER_TO_MOUSE_COUNT / $HUMAN_LIVER_COUNT * 100" | bc)
HUMAN_PANCREAS_LIFTOVER_RATE=$(echo "scale=2; $HUMAN_PANCREAS_TO_MOUSE_COUNT / $HUMAN_PANCREAS_COUNT * 100" | bc)

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE=$(echo "scale=2; $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT / $HUMAN_LIVER_TO_MOUSE_COUNT * 100" | bc)
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(echo "scale=2; $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT / $HUMAN_LIVER_TO_MOUSE_COUNT * 100" | bc)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(echo "scale=2; $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT / $HUMAN_PANCREAS_TO_MOUSE_COUNT * 100" | bc)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE=$(echo "scale=2; $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT / $HUMAN_PANCREAS_TO_MOUSE_COUNT * 100" | bc)

# Create cross-species conservation summary
CROSS_SPECIES_SUMMARY_FILE="$CROSS_SPECIES_DIR/cross_species_conservation_summary.txt"

echo "===== CROSS-SPECIES CONSERVATION SUMMARY =====" > "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "HUMAN TO MOUSE CONSERVATION STATISTICS" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total human liver regions: $HUMAN_LIVER_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total human pancreas regions: $HUMAN_PANCREAS_COUNT" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions mapped to mouse: $HUMAN_LIVER_TO_MOUSE_COUNT ($HUMAN_LIVER_LIFTOVER_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions mapped to mouse: $HUMAN_PANCREAS_TO_MOUSE_COUNT ($HUMAN_PANCREAS_LIFTOVER_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "CONSERVATION RATES" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in mouse liver: $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT ($HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in mouse pancreas: $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT ($HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse pancreas: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT ($HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse liver: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT ($HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%)" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "CROSS-SPECIES CONSERVATION ANALYSIS" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in same tissue (mouse liver): $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions conserved in different tissue (mouse pancreas): $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in same tissue (mouse pancreas): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas regions conserved in different tissue (mouse liver): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

# ===== STEP 3: FINAL COMPARATIVE ANALYSIS =====
echo "===== COMPARING TISSUE SHARING WITH SPECIES CONSERVATION ====="

# Create the comprehensive comparison file
COMPARISON_FILE="$RESULTS_DIR/tissue_vs_species_comparison.txt"

echo "===== TISSUE SHARING VS SPECIES CONSERVATION =====" > "$COMPARISON_FILE"
echo "" >> "$COMPARISON_FILE"
echo "PROJECT QUESTION 3a:" >> "$COMPARISON_FILE"
echo "For each species, tissue combination, are more open chromatin regions open in the other tissue or the other species?" >> "$COMPARISON_FILE"
echo "" >> "$COMPARISON_FILE"

# Get the within-species tissue sharing percentages
echo "WITHIN-SPECIES TISSUE SHARING" >> "$COMPARISON_FILE"
echo "-----------------------------" >> "$COMPARISON_FILE"
echo "Human liver sharing with human pancreas: $HUMAN_LIVER_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "Human pancreas sharing with human liver: $HUMAN_PANCREAS_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "Mouse liver sharing with mouse pancreas: $MOUSE_LIVER_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "Mouse pancreas sharing with mouse liver: $MOUSE_PANCREAS_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "" >> "$COMPARISON_FILE"

# Add the cross-species conservation percentages
echo "CROSS-SPECIES CONSERVATION" >> "$COMPARISON_FILE"
echo "-----------------------------" >> "$COMPARISON_FILE"
echo "Human liver conserved in mouse liver (same tissue): $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$COMPARISON_FILE"
echo "Human liver conserved in mouse pancreas (different tissue): $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$COMPARISON_FILE"
echo "Human pancreas conserved in mouse pancreas (same tissue): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$COMPARISON_FILE"
echo "Human pancreas conserved in mouse liver (different tissue): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$COMPARISON_FILE"
echo "" >> "$COMPARISON_FILE"

# Make direct comparison to answer the question
echo "DIRECT COMPARISON" >> "$COMPARISON_FILE"
echo "-----------------------------" >> "$COMPARISON_FILE"
echo "For human liver:" >> "$COMPARISON_FILE"
echo "  - Regions shared with other tissue (human pancreas): $HUMAN_LIVER_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "  - Regions conserved in same tissue other species (mouse liver): $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$COMPARISON_FILE"
echo "  => " >> "$COMPARISON_FILE"
if (( $(echo "$HUMAN_LIVER_SHARED_PCT > $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE" | bc -l) )); then
    echo "  More regions are shared with other tissue (within-species)" >> "$COMPARISON_FILE"
else
    echo "  More regions are conserved in other species (same tissue)" >> "$COMPARISON_FILE"
fi
echo "" >> "$COMPARISON_FILE"

echo "For human pancreas:" >> "$COMPARISON_FILE"
echo "  - Regions shared with other tissue (human liver): $HUMAN_PANCREAS_SHARED_PCT%" >> "$COMPARISON_FILE"
echo "  - Regions conserved in same tissue other species (mouse pancreas): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$COMPARISON_FILE"
echo "  => " >> "$COMPARISON_FILE"
if (( $(echo "$HUMAN_PANCREAS_SHARED_PCT > $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE" | bc -l) )); then
    echo "  More regions are shared with other tissue (within-species)" >> "$COMPARISON_FILE"
else
    echo "  More regions are conserved in other species (same tissue)" >> "$COMPARISON_FILE"
fi
echo "" >> "$COMPARISON_FILE"

echo "CONCLUSION:" >> "$COMPARISON_FILE"
echo "Based on the data analysis above, we can determine if transcriptional regulatory element" >> "$COMPARISON_FILE"
echo "activity is more conserved across tissues or species. The percentages indicate whether" >> "$COMPARISON_FILE"
echo "open chromatin regions tend to maintain their state more across different tissues within" >> "$COMPARISON_FILE"
echo "the same species, or across different species in the same tissue." >> "$COMPARISON_FILE"

echo "Analysis complete. Final comparison results are in: $COMPARISON_FILE"

# Create summary of all files generated
echo "===== FILES GENERATED ====="
echo "Tissue Comparison Results:"
echo "  - $TISSUE_SUMMARY_FILE"
echo "  - $HUMAN_SHARED"
echo "  - $HUMAN_LIVER_SPECIFIC"
echo "  - $HUMAN_PANCREAS_SPECIFIC"
echo "  - $MOUSE_SHARED"
echo "  - $MOUSE_LIVER_SPECIFIC"
echo "  - $MOUSE_PANCREAS_SPECIFIC"
echo ""
echo "Cross-Species Conservation Results:"
echo "  - $CROSS_SPECIES_SUMMARY_FILE"
echo "  - $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER"
echo "  - $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS"
echo "  - $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS"
echo "  - $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER"
echo "  - $HUMAN_LIVER_NOT_CONSERVED"
echo "  - $HUMAN_PANCREAS_NOT_CONSERVED"
echo ""
echo "Final Analysis:"
echo "  - $COMPARISON_FILE"#!/bin/bash
#SBATCH --job-name=cross_species
#SBATCH --output=cross_species_%j.out
#SBATCH --error=cross_species_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G

# Load modules
module load  bedtools/2.30.0

# Create output directory
RESULTS_DIR="./results"
CROSS_SPECIES_DIR="$RESULTS_DIR/cross_species"
mkdir -p $CROSS_SPECIES_DIR

# Define input files
# Original peak files
HUMAN_LIVER_PEAKS="narrowPeak/human_liver.narrowPeak"
HUMAN_PANCREAS_PEAKS="narrowPeak/human_pancreas.narrowPeak" # Updated to pancreas
MOUSE_LIVER_PEAKS="narrowPeak/mouse_liver.narrowPeak"
MOUSE_PANCREAS_PEAKS="narrowPeak/mouse_pancreas.narrowPeak" # Updated to pancreas

# HALPER output files
HUMAN_LIVER_TO_MOUSE="halper_result/idr.optimal_peak_liver.HumanToMouse.HALPER.narrowPeak.gz"
HUMAN_PANCREAS_TO_MOUSE="halper_result/idr.optimal_peak_pancreas.HumanToMouse.HALPER.narrowPeak.gz" # Updated to pancreas

# Step 1: Convert narrowPeak files to BED format for easier handling
echo "Converting files to BED format..."

# Original peak files in BED format
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$HUMAN_LIVER_PEAKS" > "$CROSS_SPECIES_DIR/human_liver.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$HUMAN_PANCREAS_PEAKS" > "$CROSS_SPECIES_DIR/human_pancreas.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$MOUSE_LIVER_PEAKS" > "$CROSS_SPECIES_DIR/mouse_liver.bed"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' "$MOUSE_PANCREAS_PEAKS" > "$CROSS_SPECIES_DIR/mouse_pancreas.bed"

# HALPER output files in BED format
zcat "$HUMAN_LIVER_TO_MOUSE" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' > "$CROSS_SPECIES_DIR/human_liver_to_mouse.bed"
zcat "$HUMAN_PANCREAS_TO_MOUSE" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' > "$CROSS_SPECIES_DIR/human_pancreas_to_mouse.bed"

# Step 2: Identify conserved regions (human regions that map to mouse and overlap with mouse peaks)
echo "Finding conserved regions across species..."

# Human liver regions conserved in mouse liver (same tissue)
HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_liver.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_liver_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_liver.bed" -u > "$HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER"

# Human pancreas regions conserved in mouse pancreas (same tissue)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_pancreas.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_pancreas_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_pancreas.bed" -u > "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS"

# Human liver regions conserved in mouse pancreas (cross-tissue conservation)
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_pancreas.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_liver_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_pancreas.bed" -u > "$HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS"

# Human pancreas regions conserved in mouse liver (cross-tissue conservation)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_liver.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_pancreas_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_liver.bed" -u > "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER"

# Step 3: Identify non-conserved regions (human regions that map to mouse but don't overlap with mouse peaks)
echo "Finding non-conserved regions across species..."

# Human liver regions not conserved in mouse (neither liver nor pancreas)
HUMAN_LIVER_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_liver_not_conserved_in_mouse.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_liver_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_liver.bed" "$CROSS_SPECIES_DIR/mouse_pancreas.bed" -v > "$HUMAN_LIVER_NOT_CONSERVED"

# Human pancreas regions not conserved in mouse (neither liver nor pancreas)
HUMAN_PANCREAS_NOT_CONSERVED="$CROSS_SPECIES_DIR/human_pancreas_not_conserved_in_mouse.bed"
bedtools intersect -a "$CROSS_SPECIES_DIR/human_pancreas_to_mouse.bed" -b "$CROSS_SPECIES_DIR/mouse_liver.bed" "$CROSS_SPECIES_DIR/mouse_pancreas.bed" -v > "$HUMAN_PANCREAS_NOT_CONSERVED"

# Step 4: Calculate statistics
echo "Calculating cross-species conservation statistics..."

# Count regions in each category
HUMAN_LIVER_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/human_liver.bed")
HUMAN_PANCREAS_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/human_pancreas.bed")
MOUSE_LIVER_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/mouse_liver.bed")
MOUSE_PANCREAS_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/mouse_pancreas.bed")

HUMAN_LIVER_TO_MOUSE_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/human_liver_to_mouse.bed")
HUMAN_PANCREAS_TO_MOUSE_COUNT=$(wc -l < "$CROSS_SPECIES_DIR/human_pancreas_to_mouse.bed")

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT=$(wc -l < "$HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS")
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT=$(wc -l < "$HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER")

HUMAN_LIVER_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_LIVER_NOT_CONSERVED")
HUMAN_PANCREAS_NOT_CONSERVED_COUNT=$(wc -l < "$HUMAN_PANCREAS_NOT_CONSERVED")

# Calculate percentages
HUMAN_LIVER_LIFTOVER_RATE=$(echo "scale=2; $HUMAN_LIVER_TO_MOUSE_COUNT / $HUMAN_LIVER_COUNT * 100" | bc)
HUMAN_PANCREAS_LIFTOVER_RATE=$(echo "scale=2; $HUMAN_PANCREAS_TO_MOUSE_COUNT / $HUMAN_PANCREAS_COUNT * 100" | bc)

HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE=$(echo "scale=2; $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT / $HUMAN_LIVER_TO_MOUSE_COUNT * 100" | bc)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(echo "scale=2; $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT / $HUMAN_PANCREAS_TO_MOUSE_COUNT * 100" | bc)
HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE=$(echo "scale=2; $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT / $HUMAN_LIVER_TO_MOUSE_COUNT * 100" | bc)
HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE=$(echo "scale=2; $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT / $HUMAN_PANCREAS_TO_MOUSE_COUNT * 100" | bc)

# Create summary report
SUMMARY_FILE="$CROSS_SPECIES_DIR/cross_species_conservation_summary.txt"

echo "===== CROSS-SPECIES CONSERVATION SUMMARY =====" > "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "HUMAN TO MOUSE CONSERVATION STATISTICS" >> "$SUMMARY_FILE"
echo "-----------------------------" >> "$SUMMARY_FILE"
echo "Total human liver regions: $HUMAN_LIVER_COUNT" >> "$SUMMARY_FILE"
echo "Total human pancreas regions: $HUMAN_PANCREAS_COUNT" >> "$SUMMARY_FILE"
echo "Human liver regions mapped to mouse: $HUMAN_LIVER_TO_MOUSE_COUNT ($HUMAN_LIVER_LIFTOVER_RATE%)" >> "$SUMMARY_FILE"
echo "Human pancreas regions mapped to mouse: $HUMAN_PANCREAS_TO_MOUSE_COUNT ($HUMAN_PANCREAS_LIFTOVER_RATE%)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "CONSERVATION RATES" >> "$SUMMARY_FILE"
echo "-----------------------------" >> "$SUMMARY_FILE"
echo "Human liver regions conserved in mouse liver: $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_COUNT ($HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%)" >> "$SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse pancreas: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_COUNT ($HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%)" >> "$SUMMARY_FILE"
echo "Human liver regions conserved in mouse pancreas: $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_COUNT ($HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%)" >> "$SUMMARY_FILE"
echo "Human pancreas regions conserved in mouse liver: $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_COUNT ($HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "CROSS-SPECIES CONSERVATION ANALYSIS" >> "$SUMMARY_FILE"
echo "-----------------------------" >> "$SUMMARY_FILE"
echo "Human liver regions conserved in same tissue (mouse liver): $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$SUMMARY_FILE"
echo "Human liver regions conserved in different tissue (mouse pancreas): $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$SUMMARY_FILE"
echo "Human pancreas regions conserved in same tissue (mouse pancreas): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$SUMMARY_FILE"
echo "Human pancreas regions conserved in different tissue (mouse liver): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "PROJECT QUESTION 3a ANALYSIS:" >> "$SUMMARY_FILE"
echo "For each species, tissue combination, are more open chromatin regions open in the other tissue or the other species?" >> "$SUMMARY_FILE"
echo "-----------------------------" >> "$SUMMARY_FILE"
echo "Looking at human liver regions:" >> "$SUMMARY_FILE"
echo "  - Conservation in mouse liver (same tissue, different species): $HUMAN_LIVER_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$SUMMARY_FILE"
echo "  - Conservation in mouse pancreas (different tissue, different species): $HUMAN_LIVER_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Looking at human pancreas regions:" >> "$SUMMARY_FILE"
echo "  - Conservation in mouse pancreas (same tissue, different species): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_PANCREAS_RATE%" >> "$SUMMARY_FILE"
echo "  - Conservation in mouse liver (different tissue, different species): $HUMAN_PANCREAS_CONSERVED_IN_MOUSE_LIVER_RATE%" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "NOTE: To fully answer this question, you should compare these cross-species conservation rates with" >> "$SUMMARY_FILE"
echo "the within-species tissue sharing rates from your previous tissue comparison analysis." >> "$SUMMARY_FILE"

echo "Analysis complete. Results are in: $CROSS_SPECIES_DIR"
