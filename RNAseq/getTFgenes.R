#!/usr/bin/env Rscript
# Grab list of TF target genes from pathway enrichment data and plot heatmap
# NOTE need to run runManualEnrich.R beforehand

# Loads all packages in a way that allows exporting to child environments
packages <- c("openxlsx", "dplyr", "tidyr", "reshape2", "pheatmap")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR/RNAseq")

######
# DEFINE I/O
######

## INPUTS ##
rc_file <- "input/Analysis files/Keith_18July19_rnaSeq_STAR_cufflinks_gencode_merged.txt"
enrich_files <- list.files(pattern="xlsx", path="output/enrichment", full.names=TRUE)

## OUTPUTS ##
output_folder <- "output/TFgenes"

######
# PREPARE DATA
######

# Read in norm readcount data
rc <- read.delim(rc_file, h=TRUE, stringsAsFactors=FALSE)

# NOTE duplicated genes: Mia2, Hspa14, St6galnac2, Nudt8
## Deal with these properly (2020-03-19)
if (grepl("log2CPM", rc_file)) {
  rc$EntrezGene.ID <- make.unique(rc$EntrezGene.ID)
  rownames(rc) <- rc$EntrezGene.ID
  rc <- rc[,-c(1:2,18)] # remove Moffat_15_S15 column (not listed in comparisons)
} else if (grepl("cufflinks", rc_file)) {
  rc$gene_short_name <- make.unique(rc$gene_short_name)
  rownames(rc) <- rc$gene_short_name
  rc <- rc[,-c(1:4,20)]
}

# Annotate columns with sample info
colAnno <- data.frame(sample=rep(NA, ncol(rc)), row.names=colnames(rc))
colAnno[1:3, "sample"] <- "MC38-Ova intergenic"
colAnno[4:6, "sample"] <- "MC38-Ova intergenic + IFN"
colAnno[7:9, "sample"] <- "MC38-Ova Atg12-1"
colAnno[c(10,11,15), "sample"] <- "MC38-Ova Atg12-1 + IFN"
colAnno[c(13,14,12), "sample"] <- "Renca-HA intergenic"
colAnno[19:21, "sample"] <- "Renca-HA Atg12-1"
colAnno[37:39, "sample"] <- "Renca-HA intergenic 2"
colAnno[43:45, "sample"] <- "Renca-HA Fitm2-1"
colAnno[16:18, "sample"] <- "Renca-HA intergenic + IFN"
colAnno[22:24, "sample"] <- "Renca-HA Atg12-1 + IFN"
colAnno[40:42, "sample"] <- "Renca-HA intergenic + IFN 2"
colAnno[46:48, "sample"] <- "Renca-HA Fitm2-1 + IFN"
colAnno[25:27, "sample"] <- "B16 intergenic"
colAnno[31:33, "sample"] <- "B16 Fitm2-1"
colAnno[28:30, "sample"] <- "B16 intergenic + IFN"
colAnno[34:36, "sample"] <- "B16 Fitm2-1 + IFN"
colAnno[49:51, "sample"] <- "Renca-ATCC WT 12hr"
colAnno[61:63, "sample"] <- "Renca-ATCC Atg12 12hr"
colAnno[52:54, "sample"] <- "Renca-ATCC WT 48hr"
colAnno[64:66, "sample"] <- "Renca-ATCC Atg12 48hr"
colAnno[73:75, "sample"] <- "Renca-ATCC Fitm2 48hr"
colAnno[58:60, "sample"] <- "Renca-ATCC WT + TNF 12hr"
colAnno[70:72, "sample"] <- "Renca-ATCC Atg12 + TNF 12hr"
colAnno[55:57, "sample"] <- "Renca-ATCC WT + IFN 48hr"
colAnno[67:69, "sample"] <- "Renca-ATCC Atg12 + IFN 48hr"
colAnno[76:78, "sample"] <- "Renca-ATCC Fitm2-1 + IFN 48hr"

# Try another annotation
colAnno2 <- data.frame(sample=rep(NA, ncol(rc)), row.names=colnames(rc))
colAnno2[1:3, "sample"] <- "Intergenic"
colAnno2[4:6, "sample"] <- "Intergenic + IFN"
colAnno2[7:9, "sample"] <- "Atg12-1"
colAnno2[c(10,11,15), "sample"] <- "Atg12-1 + IFN"
colAnno2[c(13,14,12), "sample"] <- "Intergenic"
colAnno2[19:21, "sample"] <- "Atg12-1"
colAnno2[37:39, "sample"] <- "Intergenic"
colAnno2[43:45, "sample"] <- "Fitm2-1"
colAnno2[16:18, "sample"] <- "Intergenic + IFN"
colAnno2[22:24, "sample"] <- "Atg12-1 + IFN"
colAnno2[40:42, "sample"] <- "Intergenic + IFN"
colAnno2[46:48, "sample"] <- "Fitm2-1 + IFN"
colAnno2[25:27, "sample"] <- "Intergenic"
colAnno2[31:33, "sample"] <- "Fitm2-1"
colAnno2[28:30, "sample"] <- "Intergenic + IFN"
colAnno2[34:36, "sample"] <- "Fitm2-1 + IFN"
colAnno2[49:51, "sample"] <- "WT 12hr"
colAnno2[61:63, "sample"] <- "Atg12-1 12hr"
colAnno2[52:54, "sample"] <- "WT 48hr"
colAnno2[64:66, "sample"] <- "Atg12-1 48hr"
colAnno2[73:75, "sample"] <- "Fitm2-1 48hr"
colAnno2[58:60, "sample"] <- "WT + TNF 12hr"
colAnno2[70:72, "sample"] <- "Atg12 + TNF 12hr"
colAnno2[55:57, "sample"] <- "WT + IFN 48hr"
colAnno2[67:69, "sample"] <- "Atg12 + IFN 48hr"
colAnno2[76:78, "sample"] <- "Fitm2-1 + IFN 48hr"

# Read in enrichment excel files (upregulated enrichment set)
enrich <- lapply(enrich_files, function(x) read.xlsx(x, sheet=1))
enrich <- do.call("rbind", enrich)
res_enrich <- filter(enrich, enrichment == "enriched")

# Look for pathways relating to specific TF function
tfs <- c("Irf", "Irf1", "Irf9", "Irf7", "Nfe2I2", "Nrf2", "Nfkb")

# Initiate list to hold pathway gene info
res_list <- list()

# Run for all TFs of interest
for (tf in tfs) {
  # Search result for specific TF pathway
  res_enrich_tf <- res_enrich[grep(tf, res_enrich$pathway, ignore.case=TRUE),]
  res_enrich_tf_path <- unique(res_enrich_tf$pathway)
  # Hold results in df
  res_df <- data.frame(pathway=NA, genes=NA)
  if (length(res_enrich_tf_path)) {
    for (j in seq_along(res_enrich_tf_path)) {
      res_df[j,"pathway"] <- res_enrich_tf_path[j]
      res_enrich_df_genes <- res_enrich_tf[grep(res_enrich_tf_path[j], res_enrich_tf$pathway), "genes"]
      res_enrich_df_genes <- unique(as.character(na.omit(unique(unlist(strsplit(res_enrich_df_genes, "; "))))))
      res_df[j,"genes"] <- paste(res_enrich_df_genes, collapse=", ")
    }
    # Write out df to list
    res_list[[tf]] <- res_df
  }
}

# Write out list of TF pathway genes to file
fname_file <- sprintf("%s/table_upreg_enriched_TFgenes.xlsx", output_folder)
write.xlsx(res_list, fname_file, row.names=FALSE)

for (i in seq_along(res_list)) {
  # Keep only genes of interest
  tf_name <- names(res_list)[i]
  tf_genes <- unique(unlist(strsplit(res_list[[i]]$genes, ", ")))
  tf_genes <- paste(tf_genes, collapse="|")
  rc_tf <- rc[grep(tf_genes, rownames(rc)),]

  # Plot
  if (grepl("log2CPM", rc_file)) {
    fname_out <- sprintf("%s/plot_heatmap_%s_TFgenes_normRC_log2CPM.pdf", output_folder, tf_name)
  } else if (grepl("cufflinks", rc_file)) {
    fname_out <- sprintf("%s/plot_heatmap_%s_TFgenes_normRC_FPKM.pdf", output_folder, tf_name)
  }

  title <- sprintf("Normalized read counts of %s pathway genes across mouse samples", tf_name)
  pdf(fname_out, width=15, height=7.5)
  pheatmap(rc_tf, main=title, annotation_col=colAnno)
  pheatmap(rc_tf, main=title, annotation_col=colAnno2)
  dev.off()
}
