#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=9

module load bedtools

# Set paths
HUMAN_GTF="/ocean/projects/bio230007p/ssabata/HumanGenomeInfo/gencode.v47.annotation.gff3.gz"
MOUSE_GTF="/ocean/projects/bio230007p/ssabata/MouseGenomeInfo/gencode.vM10.annotation.gff3.gz"
OUTPUT_DIR="/ocean/projects/bio230007p/ssabata/classified_regions"

# Set peak file paths
HUMAN_LIVER_PEAKS="/ocean/projects/bio230007p/ssabata/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
HUMAN_PANCREAS_PEAKS="/ocean/projects/bio230007p/ssabata/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_LIVER_PEAKS="/ocean/projects/bio230007p/ssabata/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_PANCREAS_PEAKS="/ocean/projects/bio230007p/ssabata/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"

# Create output directories
mkdir -p $OUTPUT_DIR/human
mkdir -p $OUTPUT_DIR/mouse

# Create promoter regions (2kb upstream, 200bp downstream of TSS)
echo "Creating promoter regions..."
zcat $HUMAN_GTF | awk 'BEGIN{OFS="\t"} $3=="transcript" {
    if ($7=="+") {
        start = $4-2000;
        if (start < 0) start = 0;  # Prevent negative coordinates
        print $1, start, $4+200, "promoter", ".", $7;
    } else {
        start = $5-200;
        if (start < 0) start = 0;  # Prevent negative coordinates
        print $1, start, $5+2000, "promoter", ".", $7;
    }
}' | grep -v "^chrM" > $OUTPUT_DIR/human_promoters.bed

zcat $MOUSE_GTF | awk 'BEGIN{OFS="\t"} $3=="transcript" {
    if ($7=="+") {
        start = $4-2000;
        if (start < 0) start = 0;  # Prevent negative coordinates
        print $1, start, $4+200, "promoter", ".", $7;
    } else {
        start = $5-200;
        if (start < 0) start = 0;  # Prevent negative coordinates
        print $1, start, $5+2000, "promoter", ".", $7;
    }
}' | grep -v "^chrM" > $OUTPUT_DIR/mouse_promoters.bed

# Classify liver peaks
echo "Classifying liver peaks..."
bedtools intersect -a $HUMAN_LIVER_PEAKS -b $OUTPUT_DIR/human_promoters.bed -u > $OUTPUT_DIR/human/human_liver_promoters.bed
bedtools intersect -a $HUMAN_LIVER_PEAKS -b $OUTPUT_DIR/human_promoters.bed -v > $OUTPUT_DIR/human/human_liver_enhancers.bed

bedtools intersect -a $MOUSE_LIVER_PEAKS -b $OUTPUT_DIR/mouse_promoters.bed -u > $OUTPUT_DIR/mouse/mouse_liver_promoters.bed
bedtools intersect -a $MOUSE_LIVER_PEAKS -b $OUTPUT_DIR/mouse_promoters.bed -v > $OUTPUT_DIR/mouse/mouse_liver_enhancers.bed

# Classify pancreas peaks
echo "Classifying pancreas peaks..."
bedtools intersect -a $HUMAN_PANCREAS_PEAKS -b $OUTPUT_DIR/human_promoters.bed -u > $OUTPUT_DIR/human/human_pancreas_promoters.bed
bedtools intersect -a $HUMAN_PANCREAS_PEAKS -b $OUTPUT_DIR/human_promoters.bed -v > $OUTPUT_DIR/human/human_pancreas_enhancers.bed

bedtools intersect -a $MOUSE_PANCREAS_PEAKS -b $OUTPUT_DIR/mouse_promoters.bed -u > $OUTPUT_DIR/mouse/mouse_pancreas_promoters.bed
bedtools intersect -a $MOUSE_PANCREAS_PEAKS -b $OUTPUT_DIR/mouse_promoters.bed -v > $OUTPUT_DIR/mouse/mouse_pancreas_enhancers.bed

# Calculate statistics
echo "Calculating statistics..."
echo "Human liver promoters: $(wc -l < $OUTPUT_DIR/human/human_liver_promoters.bed)" > $OUTPUT_DIR/statistics.txt
echo "Human liver enhancers: $(wc -l < $OUTPUT_DIR/human/human_liver_enhancers.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Human pancreas promoters: $(wc -l < $OUTPUT_DIR/human/human_pancreas_promoters.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Human pancreas enhancers: $(wc -l < $OUTPUT_DIR/human/human_pancreas_enhancers.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Mouse liver promoters: $(wc -l < $OUTPUT_DIR/mouse/mouse_liver_promoters.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Mouse liver enhancers: $(wc -l < $OUTPUT_DIR/mouse/mouse_liver_enhancers.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Mouse pancreas promoters: $(wc -l < $OUTPUT_DIR/mouse/mouse_pancreas_promoters.bed)" >> $OUTPUT_DIR/statistics.txt
echo "Mouse pancreas enhancers: $(wc -l < $OUTPUT_DIR/mouse/mouse_pancreas_enhancers.bed)" >> $OUTPUT_DIR/statistics.txt
