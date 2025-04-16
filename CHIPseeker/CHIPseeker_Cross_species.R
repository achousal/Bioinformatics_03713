#############################################################
# Functional annotation of conserved and specific regions   #
# between human and mouse in liver and pancreas            #
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

# Set file paths (adjust as needed)
base_path <- "/Users/yinuoyang/Desktop/cmu/03713/results/cross_species/"

# Import BED files
# Conservation and specificity files
liver_conserved <- readPeakFile(paste0(base_path, "human_liver_conserved_in_mouse_liver.bed"))
human_liver_specific <- readPeakFile(paste0(base_path, "human_liver_not_conserved_in_mouse.bed"))
mouse_liver_specific <- readPeakFile(paste0(base_path, "mouse_liver_not_conserved_in_human.bed"))

pancreas_conserved <- readPeakFile(paste0(base_path, "human_pancreas_conserved_in_mouse_pancreas.bed"))
human_pancreas_specific <- readPeakFile(paste0(base_path, "human_pancreas_not_conserved_in_mouse.bed"))
mouse_pancreas_specific <- readPeakFile(paste0(base_path, "mouse_pancreas_not_conserved_in_human.bed"))

# Set up annotation databases
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb_mouse <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate peaks
# Use human TxDb for human regions and conserved regions (represented in human coordinates)
liver_conserved_anno <- annotatePeak(liver_conserved, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_human, annoDb="org.Hs.eg.db")
human_liver_specific_anno <- annotatePeak(human_liver_specific, tssRegion=c(-3000, 3000),
                                          TxDb=txdb_human, annoDb="org.Hs.eg.db")
pancreas_conserved_anno <- annotatePeak(pancreas_conserved, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
human_pancreas_specific_anno <- annotatePeak(human_pancreas_specific, tssRegion=c(-3000, 3000),
                                             TxDb=txdb_human, annoDb="org.Hs.eg.db")

# Use mouse TxDb for mouse-specific regions
mouse_liver_specific_anno <- annotatePeak(mouse_liver_specific, tssRegion=c(-3000, 3000),
                                          TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
mouse_pancreas_specific_anno <- annotatePeak(mouse_pancreas_specific, tssRegion=c(-3000, 3000),
                                             TxDb=txdb_mouse, annoDb="org.Mm.eg.db")

# Create a directory for results
dir.create("functional_analysis_results", showWarnings = FALSE)

# Print annotation summaries
cat("Annotation Summary for Key Regions\n")
cat("\nLiver conserved regions:\n")
print(summary(liver_conserved_anno))
cat("\nHuman liver-specific regions:\n")
print(summary(human_liver_specific_anno))
cat("\nMouse liver-specific regions:\n")
print(summary(mouse_liver_specific_anno))
cat("\nPancreas conserved regions:\n")
print(summary(pancreas_conserved_anno))
cat("\nHuman pancreas-specific regions:\n")
print(summary(human_pancreas_specific_anno))
cat("\nMouse pancreas-specific regions:\n")
print(summary(mouse_pancreas_specific_anno))

# Create annotation lists for visualization
liver_anno_list <- list(
  "Conserved" = liver_conserved_anno,
  "Human_Specific" = human_liver_specific_anno,
  "Mouse_Specific" = mouse_liver_specific_anno
)

pancreas_anno_list <- list(
  "Conserved" = pancreas_conserved_anno,
  "Human_Specific" = human_pancreas_specific_anno,
  "Mouse_Specific" = mouse_pancreas_specific_anno
)

# Visualize genomic feature distribution
# Liver genomic annotations
pdf("functional_analysis_results/liver_genomic_features_bar.pdf", width=10, height=6)
plotAnnoBar(liver_anno_list)
dev.off()

pdf("functional_analysis_results/liver_genomic_features_pie.pdf", width=15, height=5)
par(mfrow=c(1,3))
plotAnnoPie(liver_conserved_anno, main="Liver Conserved Regions")
plotAnnoPie(human_liver_specific_anno, main="Human Liver-Specific")
plotAnnoPie(mouse_liver_specific_anno, main="Mouse Liver-Specific")
dev.off()

# Pancreas genomic annotations
pdf("functional_analysis_results/pancreas_genomic_features_bar.pdf", width=10, height=6)
plotAnnoBar(pancreas_anno_list)
dev.off()

pdf("functional_analysis_results/pancreas_genomic_features_pie.pdf", width=15, height=5)
par(mfrow=c(1,3))
plotAnnoPie(pancreas_conserved_anno, main="Pancreas Conserved Regions")
plotAnnoPie(human_pancreas_specific_anno, main="Human Pancreas-Specific")
plotAnnoPie(mouse_pancreas_specific_anno, main="Mouse Pancreas-Specific")
dev.off()

# Distribution of peaks relative to TSS
# Liver
pdf("functional_analysis_results/liver_tss_distance.pdf", width=10, height=6)
plotDistToTSS(liver_anno_list, title="Liver Regions - Distance to TSS")
dev.off()

# Pancreas
pdf("functional_analysis_results/pancreas_tss_distance.pdf", width=10, height=6)
plotDistToTSS(pancreas_anno_list, title="Pancreas Regions - Distance to TSS")
dev.off()

# Convert to data frames
liver_conserved_df <- as.data.frame(liver_conserved_anno)
human_liver_specific_df <- as.data.frame(human_liver_specific_anno)
mouse_liver_specific_df <- as.data.frame(mouse_liver_specific_anno)
pancreas_conserved_df <- as.data.frame(pancreas_conserved_anno)
human_pancreas_specific_df <- as.data.frame(human_pancreas_specific_anno)
mouse_pancreas_specific_df <- as.data.frame(mouse_pancreas_specific_anno)

# Extract genes for functional analysis
liver_conserved_genes <- unique(liver_conserved_df$geneId[!is.na(liver_conserved_df$geneId)])
human_liver_specific_genes <- unique(human_liver_specific_df$geneId[!is.na(human_liver_specific_df$geneId)])
mouse_liver_specific_genes <- unique(mouse_liver_specific_df$geneId[!is.na(mouse_liver_specific_df$geneId)])
pancreas_conserved_genes <- unique(pancreas_conserved_df$geneId[!is.na(pancreas_conserved_df$geneId)])
human_pancreas_specific_genes <- unique(human_pancreas_specific_df$geneId[!is.na(human_pancreas_specific_df$geneId)])
mouse_pancreas_specific_genes <- unique(mouse_pancreas_specific_df$geneId[!is.na(mouse_pancreas_specific_df$geneId)])

# Create gene lists for GO analysis
liver_gene_list <- list(
  "Conserved" = liver_conserved_genes,
  "Human_Specific" = human_liver_specific_genes
)

pancreas_gene_list <- list(
  "Conserved" = pancreas_conserved_genes,
  "Human_Specific" = human_pancreas_specific_genes
)

mouse_liver_gene_list <- list(
  "Mouse_Liver_Specific" = mouse_liver_specific_genes
)

mouse_pancreas_gene_list <- list(
  "Mouse_Pancreas_Specific" = mouse_pancreas_specific_genes
)

# GO enrichment analysis
# Human/conserved regions (using human database)
liver_go <- compareCluster(
  liver_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

pancreas_go <- compareCluster(
  pancreas_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Mouse-specific regions (using mouse database)
mouse_liver_go <- compareCluster(
  mouse_liver_gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

mouse_pancreas_go <- compareCluster(
  mouse_pancreas_gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Save GO results
write.csv(as.data.frame(liver_go), "functional_analysis_results/liver_GO_results.csv")
write.csv(as.data.frame(pancreas_go), "functional_analysis_results/pancreas_GO_results.csv")
write.csv(as.data.frame(mouse_liver_go), "functional_analysis_results/mouse_liver_GO_results.csv")
write.csv(as.data.frame(mouse_pancreas_go), "functional_analysis_results/mouse_pancreas_GO_results.csv")

# Visualize GO results
pdf("functional_analysis_results/liver_GO_dotplot.pdf", width=12, height=10)
dotplot(liver_go, showCategory=10, title="Liver - Biological Processes")
dev.off()

pdf("functional_analysis_results/pancreas_GO_dotplot.pdf", width=12, height=10)
dotplot(pancreas_go, showCategory=10, title="Pancreas - Biological Processes")
dev.off()

pdf("functional_analysis_results/mouse_liver_GO_dotplot.pdf", width=12, height=10)
dotplot(mouse_liver_go, showCategory=10, title="Mouse Liver-Specific - Biological Processes")
dev.off()

pdf("functional_analysis_results/mouse_pancreas_GO_dotplot.pdf", width=12, height=10)
dotplot(mouse_pancreas_go, showCategory=10, title="Mouse Pancreas-Specific - Biological Processes")
dev.off()

# KEGG pathway analysis
# Human/conserved regions
liver_kegg <- compareCluster(
  liver_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

pancreas_kegg <- compareCluster(
  pancreas_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Mouse-specific regions
mouse_liver_kegg <- compareCluster(
  mouse_liver_gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

mouse_pancreas_kegg <- compareCluster(
  mouse_pancreas_gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Save KEGG results
write.csv(as.data.frame(liver_kegg), "functional_analysis_results/liver_KEGG_results.csv")
write.csv(as.data.frame(pancreas_kegg), "functional_analysis_results/pancreas_KEGG_results.csv")
write.csv(as.data.frame(mouse_liver_kegg), "functional_analysis_results/mouse_liver_KEGG_results.csv")
write.csv(as.data.frame(mouse_pancreas_kegg), "functional_analysis_results/mouse_pancreas_KEGG_results.csv")

# Visualize KEGG results
pdf("functional_analysis_results/liver_KEGG_dotplot.pdf", width=12, height=8)
dotplot(liver_kegg, showCategory=10, title="Liver - KEGG Pathways")
dev.off()

pdf("functional_analysis_results/pancreas_KEGG_dotplot.pdf", width=12, height=8)
dotplot(pancreas_kegg, showCategory=10, title="Pancreas - KEGG Pathways")
dev.off()

pdf("functional_analysis_results/mouse_liver_KEGG_dotplot.pdf", width=12, height=8)
dotplot(mouse_liver_kegg, showCategory=10, title="Mouse Liver-Specific - KEGG Pathways")
dev.off()

pdf("functional_analysis_results/mouse_pancreas_KEGG_dotplot.pdf", width=12, height=8)
dotplot(mouse_pancreas_kegg, showCategory=10, title="Mouse Pancreas-Specific - KEGG Pathways")
dev.off()

# Create a summary report with statistics and top enriched terms
sink("functional_analysis_results/functional_analysis_summary.txt")
cat("# Functional Analysis of Conserved and Specific Regions\n\n")

cat("## Peak Statistics\n")
cat("Liver conserved regions:", length(liver_conserved), "\n")
cat("Human liver-specific regions:", length(human_liver_specific), "\n")
cat("Mouse liver-specific regions:", length(mouse_liver_specific), "\n")
cat("Pancreas conserved regions:", length(pancreas_conserved), "\n")
cat("Human pancreas-specific regions:", length(human_pancreas_specific), "\n")
cat("Mouse pancreas-specific regions:", length(mouse_pancreas_specific), "\n\n")

cat("## Gene Statistics\n")
cat("Genes associated with liver conserved regions:", length(liver_conserved_genes), "\n")
cat("Genes associated with human liver-specific regions:", length(human_liver_specific_genes), "\n")
cat("Genes associated with mouse liver-specific regions:", length(mouse_liver_specific_genes), "\n")
cat("Genes associated with pancreas conserved regions:", length(pancreas_conserved_genes), "\n")
cat("Genes associated with human pancreas-specific regions:", length(human_pancreas_specific_genes), "\n")
cat("Genes associated with mouse pancreas-specific regions:", length(mouse_pancreas_specific_genes), "\n\n")

# Function to extract top terms
getTopTerms <- function(go_result, cluster_name, n=10) {
  df <- as.data.frame(go_result)
  subset_df <- df[df$Cluster == cluster_name, ]
  if(nrow(subset_df) > 0) {
    subset_df <- subset_df[order(subset_df$p.adjust), ]
    return(subset_df[1:min(n, nrow(subset_df)), c("Description", "p.adjust", "Count")])
  } else {
    return(data.frame(Description=character(0), p.adjust=numeric(0), Count=numeric(0)))
  }
}

cat("## Top Biological Processes - Liver Conserved Regions\n")
liver_conserved_terms <- getTopTerms(liver_go, "Conserved")
if(nrow(liver_conserved_terms) > 0) {
  for(i in 1:nrow(liver_conserved_terms)) {
    cat(i, ". ", liver_conserved_terms$Description[i], 
        " (p-adj = ", signif(liver_conserved_terms$p.adjust[i], 3), 
        ", genes = ", liver_conserved_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

cat("\n## Top Biological Processes - Human Liver-Specific Regions\n")
human_liver_terms <- getTopTerms(liver_go, "Human_Specific")
if(nrow(human_liver_terms) > 0) {
  for(i in 1:nrow(human_liver_terms)) {
    cat(i, ". ", human_liver_terms$Description[i], 
        " (p-adj = ", signif(human_liver_terms$p.adjust[i], 3), 
        ", genes = ", human_liver_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

cat("\n## Top Biological Processes - Mouse Liver-Specific Regions\n")
mouse_liver_terms <- getTopTerms(mouse_liver_go, "Mouse_Liver_Specific")
if(nrow(mouse_liver_terms) > 0) {
  for(i in 1:nrow(mouse_liver_terms)) {
    cat(i, ". ", mouse_liver_terms$Description[i], 
        " (p-adj = ", signif(mouse_liver_terms$p.adjust[i], 3), 
        ", genes = ", mouse_liver_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

cat("\n## Top Biological Processes - Pancreas Conserved Regions\n")
pancreas_conserved_terms <- getTopTerms(pancreas_go, "Conserved")
if(nrow(pancreas_conserved_terms) > 0) {
  for(i in 1:nrow(pancreas_conserved_terms)) {
    cat(i, ". ", pancreas_conserved_terms$Description[i], 
        " (p-adj = ", signif(pancreas_conserved_terms$p.adjust[i], 3), 
        ", genes = ", pancreas_conserved_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

cat("\n## Top Biological Processes - Human Pancreas-Specific Regions\n")
human_pancreas_terms <- getTopTerms(pancreas_go, "Human_Specific")
if(nrow(human_pancreas_terms) > 0) {
  for(i in 1:nrow(human_pancreas_terms)) {
    cat(i, ". ", human_pancreas_terms$Description[i], 
        " (p-adj = ", signif(human_pancreas_terms$p.adjust[i], 3), 
        ", genes = ", human_pancreas_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

cat("\n## Top Biological Processes - Mouse Pancreas-Specific Regions\n")
mouse_pancreas_terms <- getTopTerms(mouse_pancreas_go, "Mouse_Pancreas_Specific")
if(nrow(mouse_pancreas_terms) > 0) {
  for(i in 1:nrow(mouse_pancreas_terms)) {
    cat(i, ". ", mouse_pancreas_terms$Description[i], 
        " (p-adj = ", signif(mouse_pancreas_terms$p.adjust[i], 3), 
        ", genes = ", mouse_pancreas_terms$Count[i], ")\n", sep="")
  }
} else {
  cat("No significant enrichment found\n")
}

# Extract KEGG pathway information
cat("\n## Top KEGG Pathways\n")
for(tissue in c("Liver", "Pancreas")) {
  for(region in c("Conserved", "Human-Specific", "Mouse-Specific")) {
    cat("\n### ", tissue, " - ", region, "\n", sep="")
    
    if(tissue == "Liver" && region == "Conserved") {
      kegg_df <- as.data.frame(liver_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Conserved", ]
    } else if(tissue == "Liver" && region == "Human-Specific") {
      kegg_df <- as.data.frame(liver_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Human_Specific", ]
    } else if(tissue == "Liver" && region == "Mouse-Specific") {
      kegg_df <- as.data.frame(mouse_liver_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Mouse_Liver_Specific", ]
    } else if(tissue == "Pancreas" && region == "Conserved") {
      kegg_df <- as.data.frame(pancreas_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Conserved", ]
    } else if(tissue == "Pancreas" && region == "Human-Specific") {
      kegg_df <- as.data.frame(pancreas_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Human_Specific", ]
    } else if(tissue == "Pancreas" && region == "Mouse-Specific") {
      kegg_df <- as.data.frame(mouse_pancreas_kegg)
      kegg_subset <- kegg_df[kegg_df$Cluster == "Mouse_Pancreas_Specific", ]
    }
    
    if(nrow(kegg_subset) > 0) {
      kegg_subset <- kegg_subset[order(kegg_subset$p.adjust), ]
      for(i in 1:min(5, nrow(kegg_subset))) {
        cat("- ", kegg_subset$Description[i], 
            " (p-adj = ", signif(kegg_subset$p.adjust[i], 3), 
            ", genes = ", kegg_subset$Count[i], ")\n", sep="")
      }
    } else {
      cat("No significant enrichment found\n")
    }
  }
}

cat("\n## Comparative Analysis of Conserved vs. Species-Specific Regions\n")





# Function to safely extract annotation statistics
getAnnoStats <- function(anno_obj) {
  summ <- summary(anno_obj)
  if ("annoStat" %in% names(summ) && !is.null(summ$annoStat)) {
    return(summ$annoStat$Frequency)
  } else {
    # If the expected structure isn't found, return a message
    return("Annotation statistics not available in expected format")
  }
}

# For Liver features
cat("\n### Liver: Annotation of peaks in different genomic features\n")
cat("\nConserved Regions:\n")
print(summary(liver_conserved_anno))
cat("\nHuman-Specific Regions:\n")
print(summary(human_liver_specific_anno))
cat("\nMouse-Specific Regions:\n")
print(summary(mouse_liver_specific_anno))

# For Pancreas features
cat("\n### Pancreas: Annotation of peaks in different genomic features\n")
cat("\nConserved Regions:\n")
print(summary(pancreas_conserved_anno))
cat("\nHuman-Specific Regions:\n")
print(summary(human_pancreas_specific_anno))
cat("\nMouse-Specific Regions:\n")
print(summary(mouse_pancreas_specific_anno))











#cat("\n### Liver: Percentage of peaks in different genomic features\n")
#liver_features <- rbind(
#  Conserved = summary(liver_conserved_anno)$annoStat$Frequency,
#  Human_Specific = summary(human_liver_specific_anno)$annoStat$Frequency,
#  Mouse_Specific = summary(mouse_liver_specific_anno)$annoStat$Frequency
#)
#print(liver_features)

#cat("\n### Pancreas: Percentage of peaks in different genomic features\n")
#pancreas_features <- rbind(
#  Conserved = summary(pancreas_conserved_anno)$annoStat$Frequency,
#  Human_Specific = summary(human_pancreas_specific_anno)$annoStat$Frequency,
#  Mouse_Specific = summary(mouse_pancreas_specific_anno)$annoStat$Frequency
#)
#print(pancreas_features)
#sink()

# Print completion message
cat("Functional analysis complete! All results have been saved to 'functional_analysis_results' directory.\n")
