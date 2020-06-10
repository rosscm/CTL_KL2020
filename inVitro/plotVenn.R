#!/usr/bin/env Rscript
# Script to plot overlapping hits between Renca IFNg screen and core CTL genes

# Loads all packages in a way that allows exporting to child environments
packages <- c("openxlsx", "rio", "VennDiagram", "scales", "extrafont")
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
core_file <- "input/scores/Table_S10_coreCTLgenes.xlsx"
drug_file <- "input/scores/drugz_all.xlsx"

######
# PREPARE DATA
######

# Core CTL genes
core <- read.xlsx(core_file)
core <- core$GENE

# DrugZ (subset for MUS017_Renca_IFN_end)
drug <- import_list(drug_file)
drug_r <- drug[["MUS017_Renca_IFN_end"]]

# Vectors of supp/synth genes at fdr < 0.05
drug_r_supp <- drug_r[which(drug_r$fdr_supp < 0.05), "GENE"]
drug_r_sens <- drug_r[which(drug_r$fdr_synth < 0.05), "GENE"]

######
# PLOT
######

setwd("~/Dropbox/moffat/CTL_paper/Revised Manuscript/New data/Catherine/inVitro/drugz_venn")

cols <- c("#F6EB13", "#6D90CA", "#C0BFBF")

## SUPPRESSORS ##

cols_supp <- c("#C0BFBF", "#F6EB13")
fill_supp <- c(alpha("#C0BFBF", 0.7), alpha("#F6EB13", 0.7))

# Plot
venn.diagram(
  x = list(core, drug_r_supp),
  category.names = c("Core CTL genes" , "Suppressors"),
  filename = "plot_venn_ctl_suppressor_fdr5perc.tiff",
  output = TRUE,
  units = "px",
  imagetype = "tiff",
  height = 2000,
  width = 2000,
  resolution = 500,
  compression = "lzw",
  lwd = 1,
  col = cols_supp,
  fill = fill_supp,
  cex = 1,
  fontface = "bold",
  fontfamily = "ArialMT",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.035, 0.035),
  cat.fontfamily = "ArialMT",
  ext.text = FALSE
)

## SENSITIZERS ##

cols_sens <- c("#C0BFBF", "#6D90CA")
fill_sens <- c(alpha("#C0BFBF", 0.7), alpha("#6D90CA", 0.7))

# Plot
venn.diagram(
  x = list(core, drug_r_sens),
  category.names = c("Core CTL genes" , "Sensitizers"),
  filename = "plot_venn_ctl_sensitizer_fdr5perc.tiff",
  output = TRUE,
  units = "px",
  imagetype = "tiff",
  height = 2000,
  width = 2000,
  resolution = 500,
  compression = "lzw",
  lwd = 1,
  col = cols_sens,
  fill = fill_sens,
  cex = 1,
  fontface = "bold",
  fontfamily = "ArialMT",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.035, 0.035),
  cat.fontfamily = "ArialMT",
  ext.text = FALSE
)
