module purge
module load MEME-suite/5.4.1

FASTA=mouse_conserved_sequences.fa
MOTIF_DB=H12CORE_meme_format.meme
OUTDIR=meme_output

meme-chip -oc "$OUTDIR" \
  -db "$MOTIF_DB" \
  "$FASTA"
