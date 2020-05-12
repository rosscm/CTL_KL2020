#!/usr/bin/env Rscript
# Write enrichments file for use in Cytoscape EnrichmentMap
# NOTE need to run runManualEnrich.R beforehand

# Loads all packages in a way that allows exporting to child environments
packages <- c("rio", "dplyr")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR/RNAseq")

######
# DEFINE I/O
######

## INPUTS ##
enrich_files <- list.files(pattern="xlsx", path="output/enrichment", full.names=TRUE)
enrich_names <- gsub("enrich_|_mouse_GOBP_reactome.xlsx", "", basename(enrich_files))

## OUTPUTS ##
output_folder <- "output/em_files"

######
# WORK BEGINS
######

# Generate EM table per sample analysis
for (i in seq_along(enrich_files)) {
  res <- import_list(enrich_files[i])

  # Separate up/downregulated enrichments
  for (k in seq_along(res)) {
    res_k <- res[[k]]
    res_k_name <- names(res[k])

    # Further separate into enriched/depleted pathways
    res_k_en <- filter(res_k, enrichment == "enriched", FDR < 0.05)
    res_k_dp <- filter(res_k, enrichment == "depleted", FDR < 0.05)

    # Prepare table for EM input
    if (nrow(res_k_en)) {
      # Re-name columns
      res_k_en <- res_k_en[,c(1,1,10,11)]
      colnames(res_k_en) <- c("Geneset", "Description", "NominalP", "FDR")
      # Remove gene set annotation details from Description column
      # needed for AutoAnnotate app to properly annotate aggregated gene sets
      res_k_en$Description <- gsub("\\%.*", "", res_k_en$Description)
      # Write out
      out_file_en <- sprintf("%s/table_%s_%s_enriched.txt", output_folder, enrich_names[i], res_k_name)
      write.table(res_k_en, out_file_en, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    }

    if (nrow(res_k_dp)) {
      # Re-name columns
      res_k_dp <- res_k_dp[,c(1,1,10,11)]
      colnames(res_k_dp) <- c("Geneset", "Description", "NominalP", "FDR")
      # Remove gene set annotation details from Description column
      # needed for AutoAnnotate app to properly annotate aggregated gene sets
      res_k_dp$Description <- gsub("\\%.*", "", res_k_dp$Description)
      # Write out
      out_file_dp <- sprintf("%s/table_%s_%s_depleted.txt", output_folder, enrich_names[i], res_k_name)
      write.table(res_k_dp, out_file_dp, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    }
  }
}
