#!/usr/bin/env Rscript
# Plot toptable differential expression data & analyze mValidation gene distribution

# Loads all packages in a way that allows exporting to child environments
packages <- c("reshape2", "openxlsx", "dplyr", "fgsea", "ggplot2", "ggpubr", "ggh4x", "extrafont",  "RColorBrewer", "ggrepel")
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
top_folders <- list.dirs("input/Data_cell line")
top_files <- list.files(pattern = "topTable", path = top_folders, full.names = TRUE)
annotation_file <- "input/pathway_annotations/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt"
rc_file <- "input/Analysis files/Keith_18July19_merged_STAR_log2CPM.txt"
library_file <- "input/library_annotations/mVal_full_lib.txt"
genes_file <- "input/library_annotations/mVal_gene_list.txt"

######
# PREPARE DATA
######

## NORMALIZED READS ##

# Read in norm readcount data
rc <- read.delim(rc_file, h = TRUE, stringsAsFactors = FALSE)
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

# Read in all relevant files (toptable diff expression tables)
top <- lapply(top_files, function(x) read.delim(x, h = TRUE, stringsAsFactors = FALSE))
fnames <- gsub("topTable_|.txt", "", basename(top_files))
names(top) <- fnames

## PATHWAY INFO ##

# Pathway information (all pathways minus KEGG)
pathway <- gmtPathways(annotation_file)

## CUTOFFS ##

# Define cutoff parameters for volcano plot colouring
lfc_cut = 0.4     # logFC cutoff
pval_cut = 10e-6  # p value cutoff

######
# PLOT
######

plotVolcano <- function(top_set, pathway_set) {

  ## DEG VOLCANO ##

  # Top table comparison
  top_name <- top_set

  # Organize by cell line
  cell <- unlist(strsplit(top_name, "_"))
  if (grepl("Comp", cell[1])) {
    cell_name <- cell[2]
  } else {
    cell_name <- cell[1]
  }
  cell_name <- gsub("-", "", cell_name)

  # Generate output folder
  output_folder <- cell_name
  if (!dir.exists(output_folder)) { dir.create(output_folder) }

  # Run for single comparison at a time
  top_f <- top[[top_set]]

  # Vector of all genes in analysis
  genes <- select(top_f, EntrezGene.ID) %>%
    unlist %>%
    as.character

  # Subset for genes that meet cutoff criteria (focus on upregulated genes)
  genes_select <- top_f %>%
    filter(logFC > lfc_cut, P.Value < pval_cut) %>%
    select(EntrezGene.ID) %>%
    unlist %>%
    as.character

  # Subset for pathway of interest
  pathway_name <- grep(pathway_set, names(pathway), value = TRUE)
  pathway_genes <- pathway[[pathway_name]]

  # Label genes associated with pathway
  top_f$label <- ""
  genes_label <- intersect(genes_select, pathway_genes)
  top_f[which(top_f$EntrezGene.ID %in% genes_label), "label"] <- top_f[which(top_f$EntrezGene.ID %in% genes_label), "EntrezGene.ID"]

  # Colour points by gene regulation pattern (up/down)
  top_f$colour <- "NS"
  top_f[which(top_f$logFC > lfc_cut & top_f$P.Value < pval_cut), "colour"] <- "Upregulated"
  top_f[which(top_f$logFC < -lfc_cut & top_f$P.Value < pval_cut), "colour"] <- "Downregulated"

  # Facet levels
  top_f$colour <- factor(top_f$colour, levels = c("Upregulated", "Downregulated", "NS"))

  # Get colour palette
  cols <- brewer.pal(3, "Dark2")
  cols <- c(cols[2:3], "black")

  # Plot
  p <- ggplot(top_f, aes(x = logFC, y = -log10(P.Value))) +
        geom_point(aes(fill = colour, size = colour), shape = 21) +
        scale_fill_manual(values = cols) +
        scale_size_manual(values = c(2, 2, 0.5)) +
        geom_text_repel(aes(label = label, colour = colour), size = 1.78, fontface = "italic", segment.size = 0.1) +
        geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
        geom_vline(xintercept = lfc_cut, linetype = "dashed") +
        geom_vline(xintercept = -lfc_cut, linetype = "dashed") +
        labs(x = "Log2-foldchange", y = "-Log10P") +
        theme_bw(base_size = 5) +
        theme(text = element_text(family = "ArialMT"),
              legend.position = "bottom",
              legend.title = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())

  # Define plot names
  plot_name <- sprintf("%s/plot_volcano_%s_%s", output_folder, top_name, pathway_set)

  if (pathway_set == "GO:0034976" | pathway_set == "R-HSA-381119") {
    plot_name <- paste(plot_name, "ER_v2.pdf", sep = "_")
  }
  if (pathway_set == "GO:0006986" | pathway_set == "975138.1") {
    plot_name <- paste(plot_name, "NFKB_v2.pdf", sep = "_")
  }

  # Draw out
  pdf(plot_name, width = 5, height = 5, useDingbats = FALSE)
  print(p)
  dev.off()

  ## READS BOXPLOT ##

  # Subset rc2 data
  plot_data <- rc2[which(rc2$Gene %in% genes_label),]
  if (!nrow(plot_data)) { return(NULL) }

  # Merge with topTable logFC data for point colouring
  plot_top <- select(top_f, EntrezGene.ID, logFC)
  colnames(plot_top)[1] <- "Gene"
  plot_data <- left_join(plot_data, plot_top)

  # Define up/downregulated genes from topTable data
  plot_data[which(plot_data$logFC > lfc_cut), "logFC"] <- "Upregulated"
  plot_data[which(plot_data$logFC < -lfc_cut), "logFC"] <- "Downregulated"

  # Define different cell lines
  cond <- c("MC38-Ova", "B16",
            "Renca-HA intergenic$|Renca-HA intergenic.*.IFN$|Renca-HA.*.Atg12",
            "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Atg12",
            "Renca-HA intergenic$|Renca-HA intergenic.*.IFN$|Renca-HA.*.Fitm2",
            "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Fitm2",
            "Renca-ATCC WT.*.12hr|Renca-ATCC Atg12.*.12hr",
            "Renca-ATCC WT.*.48hr|Renca-ATCC Atg12.*.48hr",
            "Renca-ATCC WT.*.48hr|Renca-ATCC Fitm2.*.48hr")

  # Init lists to store data
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

    # Reverse plotting order for condition
    if (i == "Renca-HA intergenic 2|Renca-HA intergenic.*.IFN 2|Renca-HA.*.Atg12") {
      facets_ifn <- rev(facets_ifn)
      facets_noifn <- rev(facets_noifn)
    }

    # Set levels
    plot_data2$Analysis <- factor(plot_data2$Analysis, levels = c(facets_noifn, facets_ifn))

    # Plot
    p1 <- ggplot(plot_data2, aes(x = Sample, y = Readcounts)) +
          facet_nested(~ Analysis + Sample, scales = "free") +
          geom_jitter(aes(fill = logFC), size = 2, shape = 21, width = 0.2, alpha = 0.7) +
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
                legend.position = "none",
                legend.title = element_blank())

    # Store plots in list
    plots1[[i]] <- p1

    # New plot summarizing all samples in each condition to calculate significance
    # First inverse log2 readcounts
    plot_data3 <- plot_data2
    plot_data3$Readcounts <- 2^plot_data3$Readcounts

    # Remove sample level
    plot_data3 <- plot_data3 %>%
      group_by(Gene, Analysis, logFC) %>%
      summarise(mean_reads = mean(Readcounts, na.rm = TRUE)) %>%
      as.data.frame()

    # Separate WT vs other samples
    plot_data3_wt <- plot_data3[which(plot_data3$Analysis %in% facet_ctrl),]
    plot_data3_ko <- plot_data3[-grep(paste0(facet_ctrl, "$"), plot_data3$Analysis),]

    # Define plot facets
    facet_wt <- as.character(unique(plot_data3_wt$Analysis))
    facet_ko <- as.character(unique(plot_data3_ko$Analysis))

    # Organise dfs
    colnames(plot_data3_wt)[4] <- "mean_reads_wt"
    colnames(plot_data3_ko)[4] <- "mean_reads_ko"
    plot_data3_wt <- plot_data3_wt[,-2]

    # Merge
    plot_data4 <- left_join(plot_data3_ko, plot_data3_wt)

    # Calculate log2FC
    plot_data4$mean_reads_ko <- log2(plot_data4$mean_reads_ko)
    plot_data4$mean_reads_wt <- log2(plot_data4$mean_reads_wt)
    plot_data4$mean_reads_fc <- plot_data4$mean_reads_ko - plot_data4$mean_reads_wt

    # Set levels
    plot_data4$Analysis <- factor(plot_data4$Analysis, levels = facet_ko)

    # P val group comparisons
    sample_pairs <- combn(facet_ko, 2)
    my_comparisons <- list()
    my_comparisons[[1]] <- sample_pairs[,1]
    my_comparisons[[2]] <- sample_pairs[,2]
    my_comparisons[[3]] <- sample_pairs[,3]

    # Plot
    p2 <- ggplot(plot_data4, aes(x = Analysis, y = mean_reads_fc)) +
          geom_jitter(aes(fill = logFC), size = 2, shape = 21, width = 0.2) +
          geom_boxplot(fill = "lightgrey", alpha = 0.5, position = position_dodge(), outlier.alpha = 0) +
          theme_bw(base_size = 10) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          labs(y = "Mean log2(FC) normalized reads\n(KO +/- cytokine or WT + cytokine vs. WT)") +
          scale_x_discrete(labels = function(x) { sub("\\s", "\n", x) }) +
      #    geom_text_repel(aes(label = Gene), segment.color = "grey50", size = 2.5, force = 1) +
          stat_compare_means(label = "p.signif", method = "t.test", comparisons = my_comparisons) +
          theme(text = element_text(family = "ArialMT"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.line = element_line(colour = "black"),
                legend.position = "none",
                legend.title = element_blank())

   # Store plots in list
   plots2[[i]] <- p2

   # Output table of p values
   pval <- compare_means(mean_reads_fc ~ Analysis, plot_data4, method = "t.test")
   pvals[[i]] <- pval
  }

  # Define plot names
  plot_name <- sprintf("%s/plot_boxplot_log2CPM_%s_%s", output_folder, top_set, pathway_set)
  if (pathway_set == "GO:0034976" | pathway_set == "R-HSA-381119") {
    plot_name1 <- paste(plot_name, "ER_allSamples.pdf", sep = "_")
    plot_name2 <- paste(plot_name, "ER_meanSig.pdf", sep = "_")
  }
  if (pathway_set == "GO:0006986" | pathway_set == "975138.1") {
    plot_name1 <- paste(plot_name, "NFKB_allSamples.pdf", sep = "_")
    plot_name2 <- paste(plot_name, "NFKB_meanSig.pdf", sep = "_")
  }

  # Draw out plots1
  pdf(plot_name1, width = 8, height = 5, useDingbats = FALSE)
  print(plots1)
  dev.off()

  # Draw out plots2
  pdf(plot_name2, width = 2, height = 5, useDingbats = FALSE)
  print(plots2)
  dev.off()

  # Write out pvalues
  table_name <- sprintf("%s/table_pvalues_%s_%s", output_folder, top_set, pathway_set)

  if (pathway_set == "GO:0034976" | pathway_set == "R-HSA-381119") {
    table_name <- paste(table_name, "ER_meanSig.xlsx", sep = "_")
  }
  if (pathway_set == "GO:0006986" | pathway_set == "975138.1") {
    table_name <- paste(table_name, "NFKB_meanSig.xlsx", sep = "_")
  }

  # Write out pvalues
  pvals_out <- do.call("rbind", pvals)
  write.xlsx(pvals_out, file = table_name, row.names = FALSE)
}

# Set output directory
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/RNAseq/DEG_volcano_boxplots")

# Define pathways of interest
ER_pathway  <- "GO:0034976"
ER_pathway2 <- "R-HSA-381119"
NF_pathway  <- "GO:0006986"
NF_pathway2 <- "975138.1"

# Boxplots for certain comparisons only
## NF: Renca + TNF (topTables comp 20 + 21)
## ER: B16, Renca, RencaHA Fitm2 + IFN (topTables comp 15, 19, 24, 10)

# Run function
for (f in fnames[22:23]) {
  mapply(plotVolcano, f, c(NF_pathway, NF_pathway2))
}

for (f in fnames[c(6,27,28,37)]) {
  mapply(plotVolcano, f, c(ER_pathway, ER_pathway2))
}
