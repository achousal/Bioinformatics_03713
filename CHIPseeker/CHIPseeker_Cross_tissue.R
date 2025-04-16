#############################################################
#This script is for functional annotation of human and mouse#
#conserved regions of liver and pancreas (cross-tissue)######
#############################################################
rm(list = ls())
# Install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
# For mouse 
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene") # For mouse
BiocManager::install("org.Mm.eg.db") # Mouse gene IDs
# For human 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
# Install GO.db
BiocManager::install("GO.db")
# Install clusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
# Additional packages for visualization
BiocManager::install(c("enrichplot", "ggplot2"))

# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Import the peak data (BED file)
# Human .bed file 
human_conserve_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/human_shared_between_liver_and_pancreas.bed") ##feel free to change the file path here
human_liver_specific_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/human_liver_specific.bed")
human_pancreas_specific_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/human_pancreas_specific.bed")
# Mouse .bed file
mouse_conserve_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/mouse_shared_between_liver_and_pancreas.bed") 
mouse_liver_specific_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/mouse_liver_specific.bed")
mouse_pancreas_specific_peaks <- readPeakFile("/Users/yinuoyang/Desktop/cmu/03713/tissue_comparison/mouse_pancreas_specific.bed")

# Get a summary of peak annotations
txdb_mouse <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Human 
peakAnno_human_conserve <- annotatePeak(human_conserve_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnno_human_liver <- annotatePeak(human_liver_specific_peaks, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnno_human_pancreas <- annotatePeak(human_pancreas_specific_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
# Mouse 
peakAnno_mouse_conserve <- annotatePeak(mouse_conserve_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
peakAnno_mouse_liver <- annotatePeak(mouse_liver_specific_peaks, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
peakAnno_mouse_pancreas <- annotatePeak(mouse_pancreas_specific_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_mouse, annoDb="org.Mm.eg.db")

# View annotation summary
summary(peakAnno_human_conserve)
summary(peakAnno_human_liver)
summary(peakAnno_human_pancreas)
summary(peakAnno_mouse_conserve)
summary(peakAnno_mouse_liver)
summary(peakAnno_mouse_pancreas)

# Convert to data frame for easier manipulation
peakAnnoDF_human_conserve <- as.data.frame(peakAnno_human_conserve)
peakAnnoDF_human_liver <- as.data.frame(peakAnno_human_liver)
peakAnnoDF_human_pancreas <- as.data.frame(peakAnno_human_pancreas)
peakAnnoDF_mouse_conserve <- as.data.frame(peakAnno_mouse_conserve)
peakAnnoDF_mouse_liver <- as.data.frame(peakAnno_mouse_liver)
peakAnnoDF_mouse_pancreas <- as.data.frame(peakAnno_mouse_pancreas)

# Create a list of all annotated peaks for combined visualization
all_peaks <- list(
  "Human_Conserved" = peakAnno_human_conserve,
  "Human_Liver" = peakAnno_human_liver,
  "Human_Pancreas" = peakAnno_human_pancreas,
  "Mouse_Conserved" = peakAnno_mouse_conserve,
  "Mouse_Liver" = peakAnno_mouse_liver,
  "Mouse_Pancreas" = peakAnno_mouse_pancreas
)

# Combined pie charts of genomic feature distribution
pdf("combined_pie_charts.pdf", width=12, height=8)
par(mfrow=c(2,3))
plotAnnoPie(peakAnno_human_conserve, main="Human Conserved")
plotAnnoPie(peakAnno_human_liver, main="Human Liver-Specific")
plotAnnoPie(peakAnno_human_pancreas, main="Human Pancreas-Specific")
plotAnnoPie(peakAnno_mouse_conserve, main="Mouse Conserved")
plotAnnoPie(peakAnno_mouse_liver, main="Mouse Liver-Specific")
plotAnnoPie(peakAnno_mouse_pancreas, main="Mouse Pancreas-Specific")
dev.off()

# Combined bar plot of genomic feature distribution
pdf("combined_anno_bar.pdf", width=10, height=6)
plotAnnoBar(all_peaks)
dev.off()

# Combined distribution of peaks relative to TSS
pdf("combined_distToTSS.pdf", width=10, height=6)
plotDistToTSS(all_peaks, title="Distribution of peaks relative to TSS")
dev.off()

# Extract genes for GO analysis
# Human
genes_human_conserve <- unique(peakAnnoDF_human_conserve$geneId[!is.na(peakAnnoDF_human_conserve$geneId)])
genes_human_liver <- unique(peakAnnoDF_human_liver$geneId[!is.na(peakAnnoDF_human_liver$geneId)])
genes_human_pancreas <- unique(peakAnnoDF_human_pancreas$geneId[!is.na(peakAnnoDF_human_pancreas$geneId)])

# Mouse
genes_mouse_conserve <- unique(peakAnnoDF_mouse_conserve$geneId[!is.na(peakAnnoDF_mouse_conserve$geneId)])
genes_mouse_liver <- unique(peakAnnoDF_mouse_liver$geneId[!is.na(peakAnnoDF_mouse_liver$geneId)])
genes_mouse_pancreas <- unique(peakAnnoDF_mouse_pancreas$geneId[!is.na(peakAnnoDF_mouse_pancreas$geneId)])

# Create list for comparative analysis
gene_lists_human <- list(
  Conserved = genes_human_conserve,
  Liver = genes_human_liver,
  Pancreas = genes_human_pancreas
)

gene_lists_mouse <- list(
  Conserved = genes_mouse_conserve,
  Liver = genes_mouse_liver,
  Pancreas = genes_mouse_pancreas
)

# Comparative GO enrichment analysis for human
comp_go_human <- compareCluster(
  gene_lists_human,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Comparative GO enrichment analysis for mouse
comp_go_mouse <- compareCluster(
  gene_lists_mouse,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Save results to CSV
write.csv(as.data.frame(comp_go_human), "human_GO_enrichment_results.csv")
write.csv(as.data.frame(comp_go_mouse), "mouse_GO_enrichment_results.csv")

# Visualize comparative results
pdf("human_GO_dotplot.pdf", width=12, height=10)
dotplot(comp_go_human, showCategory=15, title="Human Biological Processes Comparison")
dev.off()

pdf("mouse_GO_dotplot.pdf", width=12, height=10)
dotplot(comp_go_mouse, showCategory=15, title="Mouse Biological Processes Comparison")
dev.off()


# KEGG pathway analysis
# Convert gene IDs for KEGG (might need additional steps depending on ID format)
# Human
kegg_human <- compareCluster(
  gene_lists_human,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Mouse
kegg_mouse <- compareCluster(
  gene_lists_mouse,
  fun = "enrichKEGG",
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Save KEGG results
write.csv(as.data.frame(kegg_human), "human_KEGG_pathway_results.csv")
write.csv(as.data.frame(kegg_mouse), "mouse_KEGG_pathway_results.csv")

# KEGG pathway visualization
pdf("human_KEGG_dotplot.pdf", width=12, height=8)
dotplot(kegg_human, showCategory=15, title="Human KEGG Pathway Comparison")
dev.off()

pdf("mouse_KEGG_dotplot.pdf", width=12, height=8)
dotplot(kegg_mouse, showCategory=15, title="Mouse KEGG Pathway Comparison")
dev.off()

# Combined cross-species analysis (optional, for advanced comparison)
# We need to map mouse genes to human orthologs first

# Compare overlapping biological processes
# Create a venn diagram of enriched terms
library(VennDiagram)
library(grid)

# Extract significant GO terms
human_terms <- unique(as.data.frame(comp_go_human)$Description)
mouse_terms <- unique(as.data.frame(comp_go_mouse)$Description)

# Create Venn diagram
pdf("human_mouse_GO_overlap.pdf", width=8, height=8)
venn.plot <- draw.pairwise.venn(
  length(human_terms), 
  length(mouse_terms), 
  length(intersect(human_terms, mouse_terms)),
  category = c("Human", "Mouse"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cat.pos = c(0, 0),
  cat.dist = c(0.025, 0.025)
)
grid.draw(venn.plot)
dev.off()

# Create summary report
sink("functional_annotation_summary.txt")
cat("# Functional Annotation Analysis Summary\n\n")
cat("## Peak Statistics\n")
cat("Human conserved regions:", length(human_conserve_peaks), "\n")
cat("Human liver-specific regions:", length(human_liver_specific_peaks), "\n")
cat("Human pancreas-specific regions:", length(human_pancreas_specific_peaks), "\n")
cat("Mouse conserved regions:", length(mouse_conserve_peaks), "\n")
cat("Mouse liver-specific regions:", length(mouse_liver_specific_peaks), "\n")
cat("Mouse pancreas-specific regions:", length(mouse_pancreas_specific_peaks), "\n\n")

cat("## Gene Statistics\n")
cat("Human conserved genes:", length(genes_human_conserve), "\n")
cat("Human liver-specific genes:", length(genes_human_liver), "\n")
cat("Human pancreas-specific genes:", length(genes_human_pancreas), "\n")
cat("Mouse conserved genes:", length(genes_mouse_conserve), "\n")
cat("Mouse liver-specific genes:", length(genes_mouse_liver), "\n")
cat("Mouse pancreas-specific genes:", length(genes_mouse_pancreas), "\n\n")

cat("## Top enriched biological processes\n")
cat("### Human conserved:\n")
print(head(as.data.frame(comp_go_human)[as.data.frame(comp_go_human)$Cluster == "Conserved", c("Description", "pvalue", "qvalue")], 10))
cat("\n### Human liver-specific:\n")
print(head(as.data.frame(comp_go_human)[as.data.frame(comp_go_human)$Cluster == "Liver", c("Description", "pvalue", "qvalue")], 10))
cat("\n### Human pancreas-specific:\n")
print(head(as.data.frame(comp_go_human)[as.data.frame(comp_go_human)$Cluster == "Pancreas", c("Description", "pvalue", "qvalue")], 10))

cat("\n### Mouse conserved:\n")
print(head(as.data.frame(comp_go_mouse)[as.data.frame(comp_go_mouse)$Cluster == "Conserved", c("Description", "pvalue", "qvalue")], 10))
cat("\n### Mouse liver-specific:\n")
print(head(as.data.frame(comp_go_mouse)[as.data.frame(comp_go_mouse)$Cluster == "Liver", c("Description", "pvalue", "qvalue")], 10))
cat("\n### Mouse pancreas-specific:\n")
print(head(as.data.frame(comp_go_mouse)[as.data.frame(comp_go_mouse)$Cluster == "Pancreas", c("Description", "pvalue", "qvalue")], 10))
sink()

# Create a heatmap of top shared processes between species
# Extract common terms between human and mouse
common_terms <- intersect(human_terms, mouse_terms)

if(length(common_terms) > 0) {
  # Filter the results to get only common terms
  human_common <- as.data.frame(comp_go_human)[as.data.frame(comp_go_human)$Description %in% common_terms, ]
  mouse_common <- as.data.frame(comp_go_mouse)[as.data.frame(comp_go_mouse)$Description %in% common_terms, ]
  
  # Create a matrix for heatmap
  common_matrix <- matrix(0, nrow=length(common_terms), ncol=6)
  rownames(common_matrix) <- common_terms
  colnames(common_matrix) <- c("Human_Conserved", "Human_Liver", "Human_Pancreas", 
                               "Mouse_Conserved", "Mouse_Liver", "Mouse_Pancreas")
  
  # Fill the matrix with -log10(pvalue)
  for(i in 1:length(common_terms)) {
    term <- common_terms[i]
    # Human
    human_term <- human_common[human_common$Description == term, ]
    for(cluster in unique(human_term$Cluster)) {
      col_idx <- which(colnames(common_matrix) == paste0("Human_", cluster))
      common_matrix[i, col_idx] <- -log10(human_term$pvalue[human_term$Cluster == cluster])
    }
    
    # Mouse
    mouse_term <- mouse_common[mouse_common$Description == term, ]
    for(cluster in unique(mouse_term$Cluster)) {
      col_idx <- which(colnames(common_matrix) == paste0("Mouse_", cluster))
      common_matrix[i, col_idx] <- -log10(mouse_term$pvalue[mouse_term$Cluster == cluster])
    }
  }
  
  # Create heatmap
  pdf("cross_species_GO_heatmap.pdf", width=12, height=15)
  heatmap(common_matrix, Rowv=NA, Colv=NA, scale="none", 
          col=colorRampPalette(c("white", "red"))(100),
          margins=c(10,30))
  dev.off()
}

# Print completion message
cat("Analysis complete! All results have been saved to files in the current directory.\n")
