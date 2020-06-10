#!/usr/bin/env Rscript
# Script to analyze mValidation gene distribution

# Loads all packages in a way that allows exporting to child environments
packages <- c("reshape2", "openxlsx", "plyr", "dplyr", "fgsea", "ggplot2", "ggrepel", "ggpubr", "extrafont")
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
fc_folder <- "input/foldchange"
fc_files <- list.files(pattern = "xlsx", path = fc_folder, full.names = TRUE)
library_file <- "input/library_annotations/mVal_full_lib.txt"
genes_file <- "input/library_annotations/mVal_gene_list.txt"
score_file <- "input/library_annotations/Core_MeanNormZ.txt"

######
# PREPARE DATA
######

# Foldchange matrices (sheet 2 = guide-level mean; 3 = gene-level mean)
dat <- lapply(fc_files, function(x) read.xlsx(x, sheet = 3))

# Get query names
dat_names <- substr(basename(fc_files), 8, nchar(basename(fc_files))-40)

# mVal library gene key + z score data (defines suppressors / Sensitizer)
full_lib <- read.delim(library_file, h = TRUE, stringsAsFactors = FALSE)
gene_list <- read.delim(genes_file, h = TRUE, stringsAsFactors = FALSE)
z_score <- read.delim(score_file, h = TRUE, stringsAsFactors = FALSE)
z_score <- z_score[,-c(1,9)]
names(z_score)[2:7] <- dat_names[c(1,4,2,3,5,6)]

# Restructure gene key list for easy merging with foldchange matrices
gene_key <- melt(gene_list, measure.vars=colnames(gene_list))
colnames(gene_key) <- c("CLASS", "GENE")

# Remove empty rows
gene_key[which(gene_key$GENE == ""),] <- NA
gene_key <- na.omit(gene_key)

# Fix mismatching gene names / symbols (matching to foldchange matrices)
gene_key$GENE[gene_key$GENE == "Luciferase"] <- "luciferase"
gene_key$GENE[gene_key$GENE == "Intergenic"] <- "intergenic"
gene_key$GENE[gene_key$GENE == "CD274"] <- "Cd274"
gene_key$GENE[gene_key$GENE == "CD47"] <- "Cd47"

# Merge with z score data (z>0 = supressor; z<0 = sensitizer)
z_melt <- melt(z_score, id.vars="GENE", measure.vars=dat_names)
z_melt$CLASS <- "Non-hit"
z_melt[which(z_melt$value > 0),]$CLASS <- "Suppressor"
z_melt[which(z_melt$value < 0),]$CLASS <- "Sensitizer"

# Merge foldchange data with gene keys (eg core, control, intergenic, etc)
dat_key <- lapply(dat, function(x) join(x, gene_key, by="GENE"))

for (i in 1:length(dat_key)) {
  dat_key[[i]]$SCREEN <- NA
  dat_key[[i]]$SCREEN <- dat_names[i]

  # further classifying core genes (suppressor vs. Sensitizer)
  dat_key[[i]]$CLASS <- as.character(dat_key[[i]]$CLASS)
  z_melt_dat <- z_melt[z_melt$variable == dat_names[i],]

  for (x in which(dat_key[[i]]$GENE %in% z_melt_dat$GENE)) {
    gene <- as.character(dat_key[[i]][x,]$GENE)
    dat_key[[i]][dat_key[[i]]$GENE==gene,]$CLASS <- z_melt_dat[z_melt_dat$GENE==gene,]$CLASS
  }

  # Calculating differential log2 fold change (dropout - control)
  ## Defining mid vs late timepoints per screen
  ## Setting up rules manually; can be automated, but will take time
  dat_cols <- data.frame(tp1=c("T15.*drop.out.", "T15_Control"),
                         tp2=c("T22.*drop.out", "T24_Control"),
                         tp3=c("T11.*drop.out", "T11_Control"),
                         tp4=c("T12.*drop.out", "T12_Control"))

  if (dat_names[i] == "MC38-OVA") {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp3"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp3"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp4"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp4"], colnames(dat_key[[i]]))]
  } else if (dat_names[i] == "B16-OVA") {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp4"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp4"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp1"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp1"], colnames(dat_key[[i]]))]
  } else {
    dat_key[[i]]$Tmid_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp1"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp1"], colnames(dat_key[[i]]))]
    dat_key[[i]]$Tlate_dlogFC <- dat_key[[i]][,grep(dat_cols[1,"tp2"], colnames(dat_key[[i]]))] - dat_key[[i]][,grep(dat_cols[2,"tp2"], colnames(dat_key[[i]]))]
  }

  # Remove control + drop out columns
  cols_del <- grep(paste("*Control","*drop.out*", sep="|"), colnames(dat_key[[i]]))
  dat_key[[i]] <- dat_key[[i]][,-cols_del]
}

# Combine all screen data
dat_all <- do.call("rbind", dat_key)
dat_all <- na.omit(dat_all)

######
# PLOT
######

plotStuff <- function(tp, ptest) {

  # Only keep suppressor, sensitizer, targeting control columns
  dat_filt <- dat_all[which(dat_all$CLASS %in% c("Suppressor", "Sensitizer", "Targeting_Controls")),]
  dat_filt[dat_filt$CLASS == "Targeting_Controls",]$CLASS <- "Targeting control"
  dat_filt$CLASS <- factor(dat_filt$CLASS, levels=c("Suppressor", "Sensitizer", "Targeting control"))

  # Rename Tmid_dlogFC and Tlate_dlogFC columns
  colnames(dat_filt)[4:5] <- c("Mid", "Late")
  dat_filt$SCREEN <- factor(dat_filt$SCREEN, levels = dat_names)

  # Melt data by time point and filter by either 'mid' or 'late'
  dat_plot <- melt(dat_filt, measure.vars = c("Mid", "Late"))
  dat_plot_tp <- filter(dat_plot, variable == tp)

  # Define boxplot colours and title based on whether plotting suppressor / sensitizer
  # vs control (CLASS) or autophagy vs control (PATHWAY)
  plot_name <- sprintf("plot_boxplot_LFC_aggregated_%s_%s_v2.pdf", tp, ptest)
  table_name <- sprintf("table_pvalues_LFC_aggregated_%s_%s.txt", tp, ptest)

  # Style from boxplot_CTLpaper.R
  p <- ggplot(dat_plot_tp, aes(x = CLASS, y = value)) +
          facet_wrap(.~SCREEN, scales = "free", nrow = 1, strip.position = "bottom") +
          geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
          geom_jitter(aes(fill = CLASS), size = 1.5, shape = 21, width = 0.2) +
          geom_boxplot(aes(fill = CLASS), alpha = 0.5, outlier.alpha = 0) +
          labs(y = expression(paste(Delta, " log"[2], "-foldchange"))) +
          theme_bw(base_size = 10) +
          scale_fill_manual(values = c("#F6EB13", "#6D90CA", "#C0BFBF")) +
          theme(text = element_text(family = "ArialMT"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line = element_line(colour = "black"),
                strip.background = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank())

  # Calculate p value
  my_comparisons <- list(c("Targeting control", "Sensitizer"), c("Targeting control", "Suppressor"))
  p <- p + stat_compare_means(label = "p.signif", method = ptest, comparisons = my_comparisons)

  # P value calculations
  pval <- compare_means(value ~ CLASS, dat_plot_tp, group.by = "SCREEN", ref.group = "Targeting control", method = ptest)
  write.table(pval, file = table_name, col = TRUE, row = FALSE, quote = FALSE, sep = "\t")

  # Draw out
  pdf(plot_name, width = 10, height = 5, useDingbats = FALSE)
  print(p)
  dev.off()
}

# Plot class boxplots (sensitizer / suppressor vs control)
setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVitro/LFC_boxplots")

mapply(plotStuff, "Mid", c("t.test", "wilcox.test"))
mapply(plotStuff, "Late", c("t.test", "wilcox.test"))
