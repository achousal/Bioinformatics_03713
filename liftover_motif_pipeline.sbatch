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

# --- Step 1: Run halLiftover ---
echo "Running halLiftover..."
halLiftover /ocean/projects/bio230007p/achousal/ikaplow/Alignments/10plusway-master.hal \
  Human \
  /ocean/projects/bio230007p/achousal/Human_Liver/idr.optimal_peak.narrowPeak \
  Mouse \
  /ocean/projects/bio230007p/achousal/mouse_peaks_liftover.bed

# --- Step 2: Intersect with mouse liver peaks ---
echo "Running bedtools intersect..."
bedtools intersect -a /ocean/projects/bio230007p/achousal/mouse_peaks_liftover.bed \
  -b /ocean/projects/bio230007p/achousal/Mouse_Liver/mouse_liver_peaks.bed \
  -u > /ocean/projects/bio230007p/achousal/mouse_peaks_liftover_intersect_mouse_liver.bed

# --- Step 3: Extract FASTA sequences ---
echo "Extracting FASTA sequences..."
bedtools getfasta -fi /ocean/projects/bio230007p/achousal/ikaplow/MouseGenomeInfo/mm10.fa \
  -bed /ocean/projects/bio230007p/achousal/mouse_peaks_liftover_intersect_mouse_liver.bed \
  -fo /ocean/projects/bio230007p/achousal/mouse_conserved_sequences.fa

# --- Step 4: Run MEME-ChIP ---
echo "Running MEME-ChIP..."
meme-chip -oc /ocean/projects/bio230007p/achousal/meme_output \
  -db /ocean/projects/bio230007p/achousal/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme \
  -meme-maxw 20 \
  -meme-nmotifs 5 \
  -meme-minw 6 \
  /ocean/projects/bio230007p/achousal/mouse_conserved_sequences.fa

echo "Pipeline complete."
