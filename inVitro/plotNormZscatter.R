#!/usr/bin/env Rscript
# Plot normZ GI sector scatterplots

# Loads all packages in a way that allows exporting to child environments
packages <- c("rio", "dplyr", "ggplot2", "ggrepel", "ggExtra", "extrafont")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Import extra fonts
#font_import()
loadfonts()

# Set working directory
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/tables")

######
# DEFINE I/O
######

## INPUTS ##
dat_file1 <- normalizePath("Table_S15_Fitm2 IFN screen_DrugZ.xlsx")
dat_file2 <- normalizePath("Table_S18_Atg12 TNF screen.xlsx")

######
# FUNCTION
######

plotScatter <- function(dat_file, tp) {

  # Read specified file
  dat <- import_list(dat_file)

  # Clean tables
  for (i in seq_along(dat)) {
    colnames(dat[[i]])[-1] <- paste(colnames(dat[[i]])[-1], names(dat[i]), sep = "_")
    colnames(dat[[i]]) <- gsub("-", "_", colnames(dat[[i]]))
  }

  # Select data for specified timepoint
  sheets <- grep(tp, names(dat))
  dat_tp <- left_join(dat[[sheets[1]]], dat[[sheets[2]]], by = "GENE")
  dat_tp <- na.omit(dat_tp)

  # Define normZ columns
  norm_cols <- grep("normZ", colnames(dat_tp), value = TRUE)
  norm_wt <- grep("WT", norm_cols, value = TRUE)
  norm_cell <- setdiff(norm_cols, norm_wt)

  # Define fdr columns
  fdr_synth_cols <- grep("fdr_synth", colnames(dat_tp), value = TRUE)
  fdr_supp_cols <- grep("fdr_supp", colnames(dat_tp), value = TRUE)
  fdr_synth_wt <- grep("WT", fdr_synth_cols, value = TRUE)
  fdr_synth_cell <- setdiff(fdr_synth_cols, fdr_synth_wt)
  fdr_supp_wt <- grep("WT", fdr_supp_cols, value = TRUE)
  fdr_supp_cell <- setdiff(fdr_supp_cols, fdr_supp_wt)

  # Calculate differential normZ
  dat_tp$diff <- dat_tp[[norm_cell]] - dat_tp[[norm_wt]]
  dat_tp <- dat_tp[order(dat_tp$diff, decreasing = TRUE),]

  # Define sector lines for plotting
  hline1 <- filter(dat_tp, get(fdr_synth_cell) < 0.05) %>% summarise(median(get(norm_cell))) %>% as.numeric
  hline2 <- filter(dat_tp, get(fdr_supp_cell) < 0.05) %>% summarise(median(get(norm_cell))) %>% as.numeric
  vline1 <- filter(dat_tp, get(fdr_synth_wt) < 0.05) %>% summarise(median(get(norm_wt))) %>% as.numeric
  vline2 <- filter(dat_tp, get(fdr_supp_wt) < 0.05) %>% summarise(median(get(norm_wt))) %>% as.numeric

  # Define sector genes
  dat_tp$group <- NA
  dat_tp[which(dat_tp[[norm_cell]] >= hline1 & dat_tp[[norm_cell]] <= hline2 & dat_tp[[norm_wt]] >= vline1 & dat_tp[[norm_wt]] <= vline2), "group"] <- "middle"
  dat_tp[which(dat_tp[[norm_cell]] > hline1 & dat_tp[[norm_cell]] < hline2 & dat_tp[[norm_wt]] < vline1), "group"] <- "middle_left"
  dat_tp[which(dat_tp[[norm_cell]] > hline1 & dat_tp[[norm_cell]] < hline2 & dat_tp[[norm_wt]] > vline2), "group"] <- "middle_right"
  dat_tp[which(dat_tp[[norm_cell]] < hline1 & dat_tp[[norm_wt]] > vline1 & dat_tp[[norm_wt]] < vline2), "group"] <- "bottom_middle"
  dat_tp[which(dat_tp[[norm_cell]] < hline1 & dat_tp[[norm_wt]] < vline1), "group"] <- "bottom_left"
  dat_tp[which(dat_tp[[norm_cell]] < hline1 & dat_tp[[norm_wt]] > vline2), "group"] <- "bottom_right"
  dat_tp[which(dat_tp[[norm_cell]] > hline2 & dat_tp[[norm_wt]] > vline1 & dat_tp[[norm_wt]] < vline2), "group"] <- "top_middle"
  dat_tp[which(dat_tp[[norm_cell]] > hline2 & dat_tp[[norm_wt]] > vline2), "group"] <- "top_right"
  dat_tp[which(dat_tp[[norm_cell]] > hline2 & dat_tp[[norm_wt]] < vline1), "group"] <- "top_left"

  # Define sector labels
  dat_tp$label <- ""
  if (grepl("S15", dat_file)) {
    gene_label <- c(
      "Atg7", "Atg5", "Atg12", "Atg10", "Jak1", "Jak2", "Atg3", "Pex19", "Fitm2",
      "Ifngr1", "Ifngr2", "Stat1", "Pex1", "Pex16", "Pex3", "Pex6", "Pex16",
      "Ptpn2", "Stub1", "Gclc", "Dhx36"
    )
    dat_tp[which(dat_tp$GENE %in% gene_label), "label"] <- dat_tp[which(dat_tp$GENE %in% gene_label), "GENE"]
  } else if (grepl("S18", dat_file)) {
    gene_label <- c(
      "Tnfrsf1a", "Tradd", "Sqstm1", "Atg16l1", "Atg5", "Atg7", "Atg12", "Atg3",
      "Atg10", "Tbk1", "Ripk1", "Atg101", "Ikbkg", "Atg9a", "Kdm8", "Atg14", "Traf2",
      "Tnfaip3", "Ikbkb", "Traf6", "Pik3c3", "Atg2a", "Vps37a", "Fbxw11", "Men1"
    )
    dat_tp[which(dat_tp$GENE %in% gene_label), "label"] <- dat_tp[which(dat_tp$GENE %in% gene_label), "GENE"]
  }

  # Define point fill
  dat_tp$type <- "none"
  dat_tp[grep("top", dat_tp$group), "type"] <- "suppressor"
  dat_tp[grep("bottom", dat_tp$group), "type"] <- "sensitizer"
  dat_tp$type <- factor(dat_tp$type)

  # Define plot axis labels
  x_lab <- paste(unlist(lapply(strsplit(norm_wt, "_"), "[", 2)), "(NormZ)")
  y_lab <- paste(unlist(lapply(strsplit(norm_cell, "_"), "[", 2)), "(NormZ)")

  # Plot
  p <- ggplot(dat_tp, aes(x = get(norm_wt), y = get(norm_cell), fill = type)) +
        geom_hline(yintercept = hline1, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_hline(yintercept = hline2, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_vline(xintercept = vline1, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_vline(xintercept = vline2, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_point(data = subset(dat_tp, group == "middle"), size = 0.1, shape = 21, colour = "grey80", fill = "grey80", stroke = 0.1) +
        geom_point(data = subset(dat_tp, group != "middle"), aes(size = abs(diff)), shape = 21, stroke = 0.1) +
        scale_fill_manual(values = c("grey80", "#00AEEF", "#FCCC00")) +
        scale_colour_manual(values = c("grey80", "#00AEEF", "#FCCC00")) +
        scale_size_continuous(range = c(0.5, 2)) +
        geom_text_repel(aes(label = label, colour = type), size = 1.78, fontface = "italic", segment.size = 0.1) +
        labs(x = x_lab, y = y_lab) +
        theme_bw(base_size = 5) +
        theme(text = element_text(family = "ArialMT"),
              legend.position = "none",
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())

  # With marginal distributions
  #p2 <- ggMarginal(p, groupFill = TRUE)

  # Draw out
  fname <- gsub(".xlsx", "", basename(dat_file))
  plot_name <- sprintf("plot_scatter_%s_%s_v2.pdf", fname, tp)
  plot_name <- gsub(" ", "", plot_name)
  pdf(plot_name, width = 3, height = 3, useDingbats = FALSE)
  print(p)
  dev.off()
}

# Plot scatterplots
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVitro/drugz_scatter")

# Run function on different data files
mapply(plotScatter, dat_file1, "End")
mapply(plotScatter, dat_file2, "End")
