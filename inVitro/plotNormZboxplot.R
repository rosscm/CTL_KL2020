#!/usr/bin/env Rscript
# Script to analyze mValidation gene distribution

# Loads all packages in a way that allows exporting to child environments
packages <- c("reshape2", "openxlsx", "plyr", "dplyr", "fgsea", "ggplot2", "ggrepel", "ggpubr", "KEGGREST", "extrafont")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Import extra fonts
font_import()
loadfonts()

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR")

######
# DEFINE I/O
######

## INPUTS ##
score_file <- "input/scores/Table_S7_DrugZ output.xlsx"
library_file <- "input/library_annotations/mVal_full_lib.txt"
genes_file <- "input/library_annotations/mVal_gene_list.txt"
annotation_file <- "input/pathway_annotations/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt"

######
# PREPARE DATA
######

## DRUGZ ##

# Read in table
score <- read.xlsx(score_file)
score <- score[,-1]

# Split into mid/late + normZ/FDR
## Mid normz
mid_score <- score[grep("GENE|Mid.NormZ", colnames(score))]
colnames(mid_score) <- gsub("\\..*", "", colnames(mid_score))
mid_score <- melt(mid_score)
colnames(mid_score) <- c("Gene", "Screen", "NormZ")
## Mid FDR
mid_fdr <- score[grep("GENE|Mid.FDR", colnames(score))]
colnames(mid_fdr) <- gsub("\\..*", "", colnames(mid_fdr))
mid_fdr <- melt(mid_fdr)
colnames(mid_fdr) <- c("Gene", "Screen", "FDR")
## Mid
mid <- join(mid_score, mid_fdr)
mid$Group <- "Mid"

## End normz
end_score <- score[grep("GENE|End.NormZ", colnames(score))]
colnames(end_score) <- gsub("\\..*", "", colnames(end_score))
end_score <- melt(end_score)
colnames(end_score) <- c("Gene", "Screen", "NormZ")
## End FDR
end_fdr <- score[grep("GENE|End.FDR", colnames(score))]
colnames(end_fdr) <- gsub("\\..*", "", colnames(end_fdr))
end_fdr <- melt(end_fdr)
colnames(end_fdr) <- c("Gene", "Screen", "FDR")
## End
end <- join(end_score, end_fdr)
end$Group <- "End"

## All
all_scores <- rbind(mid, end)

## LIBRARY INFO ##

# mVal library gene key
genes <- read.delim(genes_file, h=TRUE, stringsAsFactors=FALSE)

# Restructure gene key list for easy merging with foldchange matrices
gene_key <- melt(genes, measure.vars=colnames(genes))
colnames(gene_key) <- c("Type", "Gene")

# Remove empty rows
gene_key[which(gene_key$Gene == ""),] <- NA
gene_key <- na.omit(gene_key)

# Fix mismatching gene names / symbols (matching to score matrices)
gene_key$Gene[gene_key$Gene == "Luciferase"] <- "luciferase"
gene_key$Gene[gene_key$Gene == "Intergenic"] <- "intergenic"
gene_key$Gene[gene_key$Gene == "CD274"] <- "Cd274"
gene_key$Gene[gene_key$Gene == "CD47"] <- "Cd47"

## COMBINE ##

# Merge with normz data (z>0 = supressor; z<0 = sensitizer)
z_melt <- left_join(all_scores, gene_key)

# Define supp/sens
z_melt$Class <- NA
z_melt[which(z_melt$NormZ > 0),]$Class <- "Resistors"
z_melt[which(z_melt$NormZ < 0),]$Class <- "Sensitizers"
z_melt[which(z_melt$Type == "Targeting_Controls"),]$Class <- "Targeting_control"
z_melt[which(z_melt$Type == "Non.targeting.controls"),]$Class <- "Non_targeting_control"
z_melt[which(z_melt$Type == "Others"),]$Class <- "Other"

# Remove type column
z_melt <- z_melt[,-6]

# Label gene if significant (FDR < 0.05)
z_melt$Label <- ""
z_melt[which(z_melt$FDR < 0.05), "Label"] <- z_melt[which(z_melt$FDR < 0.05), "Gene"]
z_melt <- na.omit(z_melt)

## PATHWAY INFO ##

# Pathway information (all pathways minus KEGG)
pathway <- gmtPathways(annotation_file)

# KEGG ferroptosis pathway using KEGGREST
kegg <- keggGet("mmu04216")
ferr_pathway <- kegg[[1]]$ENTRY
ferr_genes <- kegg[[1]]$GENE
ferr_genes <- ferr_genes[grep(";", ferr_genes)]
ferr_genes <- unlist(lapply(strsplit(ferr_genes, ";"), "[", 1))

# Create pathway
ferr <- list()
ferr[[ferr_pathway]] <- ferr_genes

# Add to pathway list
pathway <- c(pathway, ferr)

######
# PLOT
######

plotBoxplot <- function(pathway_set) {

  # Subset for pathway of interest
  pathway_name <- grep(pathway_set, names(pathway), value = TRUE)
  pathway_genes <- pathway[[pathway_name]]

  # Subset z_melt data
  plot_data <- z_melt[which(z_melt$Gene %in% pathway_genes),]
  plot_data$Group <- factor(plot_data$Group, levels = c("Mid", "End"))

  # Only colour points if significant
  plot_data[which(plot_data$Label == ""), "Class"] <- "NS"
  plot_data$Class <- factor(plot_data$Class, levels = c("Resistors", "Sensitizers", "NS"))

  # Remove Atg genes
  plot_data <- plot_data[-grep("Atg", plot_data$Gene),]

  # Draw out
  plot_name <- sprintf("plot_boxplot_genomewide_NormZ_%s", pathway_set)

  if (pathway_set == "R-HSA-381119"| pathway_set == "GO:0034976") { plot_name <- paste(plot_name, "ER.pdf", sep = "_") }
  if (pathway_set == "mmu04216") { plot_name <- paste(plot_name, "Ferr_noAtg.pdf", sep = "_") }

  pdf(plot_name, width = 6.5, height = 4, useDingbats = FALSE)

  p <- ggplot(plot_data, aes(x = Group, y = NormZ, label = Label)) +
        facet_grid(. ~ Screen) +
        geom_jitter(aes(fill = Class), size = 1.5, shape = 21, width = 0.2) +
        geom_boxplot(fill = "lightgrey", alpha = 0.5, position = position_dodge(), outlier.alpha = 0) +
        geom_text_repel(aes(label = Label), segment.color = "grey50", size = 2.5, force = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
        theme_bw(base_size = 10) +
        scale_fill_manual(values = c("#F6EB13", "#6D90CA", "#C0BFBF")) +
        theme(text = element_text(family = "ArialMT"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
  dev.off()
}

setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New data/Catherine/inVitro/pathway_boxplots")

# Run function on pathway sets
ER_pathway <- "R-HSA-381119"
ER_pathway2 <- "GO:0034976"
Ferr_pathway <- "mmu04216"

mapply(plotBoxplot, c(ER_pathway, ER_pathway2, Ferr_pathway))
