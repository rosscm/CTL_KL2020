#!/usr/bin/env Rscript
# Script to analyze mValidation gene distribution

# Loads all packages in a way that allows exporting to child environments
packages <- c("reshape2", "openxlsx", "plyr", "dplyr", "fgsea", "ggplot2", "ggrepel", "ggpubr", "ggh4x", "extrafont")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Import extra fonts
#font_import()
loadfonts()

# Set working directory
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine")

######
# DEFINE I/O
######

## INPUTS ##
#rc_file <- "input/Analysis files/Keith_18July19_merged_STAR_readCounts.txt"
rc_file <- "input/Analysis files/Keith_18July19_merged_STAR_log2CPM.txt"
deg_folders <- list.dirs("input/Data_cell line")
deg_files <- list.files(pattern="topTable", path = deg_folders, full.names = TRUE)
library_file <- "input/library_annotations/mVal_full_lib.txt"
genes_file <- "input/library_annotations/mVal_gene_list.txt"
annotation_file <- "input/pathway_annotations/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt"

######
# PREPARE DATA
######

## NORMALIZED READS ##

# Read in norm readcount data
rc <- read.delim(rc_file, h=TRUE, stringsAsFactors=FALSE)
rc <- rc[ !apply( rc[,-1:-2] == 0,1, all), ]
genes <- rc$EntrezGene.ID
rc <- rc[,-c(1:2,18)] # remove Moffat_15_S15 column (not listed in comparisons)

# Annotate columns with sample info
colAnno <- data.frame(Sample = colnames(rc), Analysis = rep(NA, ncol(rc)))
colAnno[1:3, "Analysis"] <- "MC38-Ova intergenic"
colAnno[4:6, "Analysis"] <- "MC38-Ova intergenic + IFN"
colAnno[7:9, "Analysis"] <- "MC38-Ova Atg12-1"
colAnno[c(10,11,15), "Analysis"] <- "MC38-Ova Atg12-1 + IFN"
colAnno[c(13,14,12), "Analysis"] <- "Renca-HA intergenic"
colAnno[19:21, "Analysis"] <- "Renca-HA Atg12-1"
colAnno[37:39, "Analysis"] <- "Renca-HA intergenic 2"
colAnno[43:45, "Analysis"] <- "Renca-HA Fitm2-1"
colAnno[16:18, "Analysis"] <- "Renca-HA intergenic + IFN"
colAnno[22:24, "Analysis"] <- "Renca-HA Atg12-1 + IFN"
colAnno[40:42, "Analysis"] <- "Renca-HA intergenic + IFN 2"
colAnno[46:48, "Analysis"] <- "Renca-HA Fitm2-1 + IFN"
colAnno[25:27, "Analysis"] <- "B16 intergenic"
colAnno[31:33, "Analysis"] <- "B16 Fitm2-1"
colAnno[28:30, "Analysis"] <- "B16 intergenic + IFN"
colAnno[34:36, "Analysis"] <- "B16 Fitm2-1 + IFN"
colAnno[49:51, "Analysis"] <- "Renca-ATCC WT 12hr"
colAnno[61:63, "Analysis"] <- "Renca-ATCC Atg12 12hr"
colAnno[52:54, "Analysis"] <- "Renca-ATCC WT 48hr"
colAnno[64:66, "Analysis"] <- "Renca-ATCC Atg12 48hr"
colAnno[73:75, "Analysis"] <- "Renca-ATCC Fitm2 48hr"
colAnno[58:60, "Analysis"] <- "Renca-ATCC WT + TNF 12hr"
colAnno[70:72, "Analysis"] <- "Renca-ATCC Atg12 + TNF 12hr"
colAnno[55:57, "Analysis"] <- "Renca-ATCC WT + IFN 48hr"
colAnno[67:69, "Analysis"] <- "Renca-ATCC Atg12 + IFN 48hr"
colAnno[76:78, "Analysis"] <- "Renca-ATCC Fitm2-1 + IFN 48hr"

# Add gene column back to rc
rc <- cbind(Gene = genes, rc)

# Melt data
rc2 <- melt(rc)
colnames(rc2) <- c("Gene", "Sample", "Readcounts")

# Join with colAnno
rc2 <- left_join(rc2, colAnno)
rc2 <- rc2[,c(1,2,4,3)]

## DEG INFO ##

# Get files outlining TNF effect in WT Renca cells
## Comp 20 + 21 - simple pairwise models for Tx vs UnTx (KB)



## PATHWAY INFO ##

# Pathway information (all pathways minus KEGG)
pathway <- gmtPathways(annotation_file)

######
# PLOT
######

plotBoxplot <- function(pathway_set) {

  # Subset for pathway of interest
  pathway_name <- grep(pathway_set, names(pathway), value = TRUE)
  pathway_genes <- pathway[[pathway_name]]

  # Subset rc_melt data
  plot_data <- rc2[which(rc2$Gene %in% pathway_genes),]

  # Define different cell lines
  cond <- c("MC38-Ova", "B16",
            "Renca-HA intergenic$|Renca-HA intergenic.*.IFN$|Renca-HA.*.Atg12",
            "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Atg12",
            "Renca-HA intergenic$|Renca-HA intergenic.*.IFN$|Renca-HA.*.Fitm2",
            "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Fitm2",
            "Renca-ATCC WT.*.12hr|Renca-ATCC Atg12.*.12hr",
            "Renca-ATCC WT.*.48hr|Renca-ATCC Atg12.*.48hr",
            "Renca-ATCC WT.*.48hr|Renca-ATCC Fitm2.*.48hr")

  plots1 <- list()
  plots2 <- list()
  pvals <- list()

  for (i in cond) {

    # Sort data for plotting
    plot_data2 <- plot_data[grep(i, plot_data$Analysis),]
    facets <- unique(plot_data2$Analysis)
    facets_ifn <- grep("IFN|TNF", facets, value = TRUE)
    facets_noifn <- facets[!(facets %in% facets_ifn)]
    facet_ctrl <- grep("intergenic|WT", facets_noifn, value = TRUE)

    if (i == "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Atg12") {
      facets_ifn <- rev(facets_ifn)
      facets_noifn <- rev(facets_noifn)
    }

    # Set levels
    plot_data2$Analysis <- factor(plot_data2$Analysis, levels = c(facets_noifn, facets_ifn))

    # Plot
    p1 <- ggplot(plot_data2, aes(x = Sample, y = Readcounts)) +
          facet_nested(~ Analysis + Sample, scales = "free") +
          geom_jitter(fill = "black", size = 1.5, shape = 21, width = 0.2) +
          geom_boxplot(fill = "lightgrey", alpha = 0.5, position = position_dodge(), outlier.alpha = 0) +
          theme_bw(base_size = 10) +
          labs(y = "Log2(CPM) normalized reads") +
          theme(text = element_text(family = "sans"),
                panel.grid = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.border = element_blank(),
                axis.title.x = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position = "bottom",
                legend.title = element_blank())

    # Store plots in list
    plots1[[i]] <- p1

    # New plot summarizing all samples in each condition to calculate significance
    # Remove sample level
    plot_data3 <- plot_data2 %>%
      group_by(Gene, Analysis) %>%
      summarise(mean_reads = mean(Readcounts)) %>%
      as.data.frame()

    # Set levels
    plot_data3$Analysis <- factor(plot_data3$Analysis, levels = c(facets_noifn, facets_ifn))

    # Plot
    p2 <- ggplot(plot_data3, aes(x = Analysis, y = mean_reads)) +
          geom_jitter(fill = "black", size = 1.5, shape = 21, width = 0.2) +
          geom_boxplot(fill = "lightgrey", alpha = 0.5, position = position_dodge(), outlier.alpha = 0) +
          theme_bw(base_size = 10) +
          labs(y = "Mean log2(CPM) normalized reads") +
          scale_x_discrete(labels = function(x) { sub("\\s", "\n", x) }) +
          stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = facet_ctrl,
                             label.y = max(plot_data3$mean_reads) + 3) +
          theme(text = element_text(family = "ArialMT"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.line = element_line(colour = "black"),
                legend.position = "bottom",
                legend.title = element_blank())

   # Store plots in list
   plots2[[i]] <- p2

   # Output table of p values
   pval <- compare_means(mean_reads ~ Analysis, plot_data3, ref.group = facet_ctrl, method = "wilcox.test")
   pvals[[i]] <- pval
  }

  # Define plot names
  plot_name <- sprintf("plot_boxplot_genomewide_log2CPM_%s", pathway_set)

  if (pathway_set == "GO:0034976" | pathway_set == "R-HSA-381119") {
    plot_name1 <- paste(plot_name, "ER_allSamples.pdf", sep = "_")
    plot_name2 <- paste(plot_name, "ER_meanSig.pdf", sep = "_")
  }
  if (pathway_set == "GO:0006986" | pathway_set == "975138.1") {
    plot_name1 <- paste(plot_name, "NFKB_allSamples.pdf", sep = "_")
    plot_name2 <- paste(plot_name, "NFKB_meanSig.pdf", sep = "_")
  }

  # Draw out
  pdf(plot_name1, width = 8, height = 5, useDingbats = FALSE)
  print(plots1)
  dev.off()

  pdf(plot_name2, width = 2.5, height = 5, useDingbats = FALSE)
  print(plots2)
  dev.off()

  # Write out pvalues
  table_name <- sprintf("table_pvalues_%s", pathway_set)

  if (pathway_set == "GO:0034976" | pathway_set == "R-HSA-381119") {
    table_name <- paste(table_name, "ER_meanSig.xlsx", sep = "_")
  }
  if (pathway_set == "GO:0006986" | pathway_set == "975138.1") {
    table_name <- paste(table_name, "NFKB_meanSig.xlsx", sep = "_")
  }

  pvals_out <- do.call("rbind", pvals)
  write.xlsx(pvals_out, file = table_name, row.names = FALSE)

}

setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New data/Catherine/RNAseq/pathway_boxplots")

# Run function on pathway sets
ER_pathway  <- "GO:0034976"
ER_pathway2 <- "R-HSA-381119"
NF_pathway  <- "GO:0006986"
NF_pathway2 <- "975138.1"

mapply(plotBoxplot, c(ER_pathway, ER_pathway2, NF_pathway, NF_pathway2))
