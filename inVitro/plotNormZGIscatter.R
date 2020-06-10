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
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/input/tables")

######
# DEFINE I/O
######

## INPUTS ##
dat_file1 <- normalizePath("Table_S17_Fitm2_GI_data.xlsx")
dat_file2 <- normalizePath("Table_S19_Atg12_GI_data.xlsx")

######
# FUNCTION
######

plotScatter <- function(dat_file, tp, treatment) {

  # Read specified file
  dat <- import_list(dat_file)

  # Clean tables
  dat <- lapply(dat, function(x) {
    colnames(x)[1] <- "GENE"
    return(x)
  })

  # Select data for specified treatment + timepoint
  sheet_name <- paste(c(treatment, tp), collapse = "_")
  sheet_num <- grep(sheet_name, names(dat))
  dat_tp <- dat[[sheet_num]]
  dat_tp <- na.omit(dat_tp)

  # Define normZ columns of interest
  norm_wt <- colnames(dat_tp[2])
  norm_cell <- colnames(dat_tp[3])

  # Define plot axis labels
  x_lab <- paste(norm_wt, "(NormZ)")
  y_lab <- paste(norm_cell, "(NormZ)")

  # Define sector lines for plotting (median normZ of sens/supp in WT and mutant background)
  hline1 <- filter(dat_tp, get(norm_cell) < 0, fdr < 0.05) %>% summarise(median(get(norm_cell))) %>% as.numeric
  hline2 <- filter(dat_tp, get(norm_cell) > 0, fdr < 0.05) %>% summarise(median(get(norm_cell))) %>% as.numeric
  vline1 <- filter(dat_tp, get(norm_wt) < 0, fdr < 0.05) %>% summarise(median(get(norm_wt))) %>% as.numeric
  vline2 <- filter(dat_tp, get(norm_wt) > 0, fdr < 0.05) %>% summarise(median(get(norm_wt))) %>% as.numeric

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

  # Define point fill
  dat_tp$type <- "none"
  dat_tp[which(dat_tp[[norm_cell]] <= hline1 & dat_tp$fdr < 0.05), "type"] <- "sensitizer"
  dat_tp[which(dat_tp[[norm_cell]] >= hline2 & dat_tp$fdr < 0.05), "type"] <- "suppressor"
  dat_tp$type <- factor(dat_tp$type)

  # Define sector labels
  dat_tp$label <- ""
  if (grepl("Fitm2", dat_file) & grepl("IFN", sheet_name)) {
    gene_label <- c(
      "Atg7", "Atg12", "Atg5", "Atg16l1", "Atg10", "Atg3", "Fitm2", "Irf1",
      "Irgm1", "Irgm2", "Pex1", "Pex16", "Pex6", "Cflar", "Vps13a", "Stub1",
      "Rnf114", "Elovl1", "Parp10", "Trim32", "Pex2", "Pex3", "Pex19"
    )
  } else if (grepl("Fitm2", dat_file) & grepl("CTL", sheet_name)) {
    gene_label <- c(
      "Jak2", "Paics", "Ldha", "Inip", "Vps29", "Irf1", "Vps18", "Rnf31", "Vps16",
      "Fitm2", "Glrx3", "Hexim1", "Tap2", "B2m", "Ube2g2", "Slc39a6", "Sephs1",
      "Ept1", "Tmem41b", "Far1", "Alg8", "Atp71", "Alg5", "Trim32", "Parp10"
    )
  } else if (grepl("Atg12", dat_file) & grepl("TNF", sheet_name)) {
    gene_label <- c(
      "Tnfrsf1a", "Sqstm1", "Tradd", "Atg16l1", "Nfkbia", "Atg5", "Casp8", "Atg12",
      "Atg7", "Atg3", "Atg10", "Ripk1", "Rnf31", "Tbk1", "Ei24", "Ikbkg", "Atg101",
      "Chuk", "Men1", "Dot1", "Rab1b", "Gnas", "Kmt2a", "Gdi2", "F8a"
    )
  } else if (grepl("Atg12", dat_file) & grepl("CTL", sheet_name)) {
    gene_label <- c(
      "Rela", "B2m", "Rab1b", "Peli1", "Gnas", "Klhl9", "Tap2", "Raver1", "Men1",
      "Echs1", "Chuk", "Pcnt", "Smad4", "Rabif", "Irf2", "Ankrd49", "Fas", "Tatdn2",
      "Fbxw11", "Rbm15", "Wdr45", "Syvn1", "Ube2g2", "Pten", "Det1", "Eif2b5",
      "Jak2", "Dnaaf5", "Rheb", "Sqstm1", "Med23", "Med25", "Mta2", "Med1", "l7Rn6",
      "Snrpd3", "Tbk1", "Amd1", "Ccnc", "Atg10", "H2−T23", "Zbtb7a", "Atg7", "Med13",
      "Iigp1", "Irf9", "Atg12", "Hira", "Med12", "Irf1", "Atg5", "Atg16"
    )
  }

  # Label specified points + all coloured points
  dat_tp[which(dat_tp$GENE %in% gene_label), "label"] <- dat_tp[which(dat_tp$GENE %in% gene_label), "GENE"]
  dat_tp[which(dat_tp$type == "sensitizer" | dat_tp$type == "suppressor"), "label"] <- dat_tp[which(dat_tp$type == "sensitizer" | dat_tp$type == "suppressor"), "GENE"]

  # Override labels for genes with FDR >= 0.05
  #dat_tp[which(dat_tp[[norm_cell]] >= hline1 & dat_tp[[norm_cell]] <= hline2 & dat_tp[[norm_wt]] >= vline1 & dat_tp[[norm_wt]] <= vline2), "label"] <- ""
  dat_tp[which(dat_tp$fdr >= 0.05), "label"] <- ""

  # Plot
  p <- ggplot(dat_tp, aes(x = get(norm_wt), y = get(norm_cell), fill = type)) +
        geom_hline(yintercept = hline1, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_hline(yintercept = hline2, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_vline(xintercept = vline1, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_vline(xintercept = vline2, linetype = "dashed", size = 0.1, colour = "grey80") +
        geom_point(data = subset(dat_tp, fdr >= 0.05), size = 0.1, shape = 21, colour = "grey80", fill = "grey80", stroke = 0.1) +
        geom_point(data = subset(dat_tp, fdr < 0.05), aes(size = abs(log2(fdr))), shape = 21, stroke = 0.1) +
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
  plot_name <- sprintf("plot_scatter_%s_%s_%s.pdf", fname, tp, treatment)
  plot_name <- gsub(" ", "", plot_name)
  pdf(plot_name, width = 3, height = 3, useDingbats = FALSE)
  print(p)
  dev.off()
}

# Set output directory
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVitro/drugz_scatter")

# Run function on different data files
mapply(plotScatter, dat_file1, c("en", "mi"), "IFN")
mapply(plotScatter, dat_file1, c("en", "mi"), "CTL")
mapply(plotScatter, dat_file2, c("en", "mi"), "TNF")
mapply(plotScatter, dat_file2, c("en", "mi"), "CTL")
