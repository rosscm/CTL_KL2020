#!/usr/bin/env Rscript
# Run GO enrichment on various sets differentially expressed genes

# Loads all packages in a way that allows exporting to child environments
packages <- c("fgsea", "plyr", "dplyr", "openxlsx")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR/RNAseq")

######
# DEFINE I/O
######

## INPUTS ##
main_folder <- "input/Data_cell line"
data_folders <- list.dirs(main_folder)
annotation_file <- "input/pathway_annotations/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt"

## INPUTS ##
output_folder <- "output/enrichment"

######
# PARAMETER SETTING
######

path_min = 10   # min pathway size
path_max = 300  # max pathway size
lfc_cut = 0.4   # logFC cutoff
fdr_cut = 0.05  # FDR cutoff
exp_cut = 0     # average expression cutoff

######
# PREPARE DATA
######

# List all relevant files (toptable diff expression tables)
top_files <- list.files(pattern=".txt", path=data_folders, full.names=TRUE)
top <- lapply(top_files, function(x) read.delim(x, h=TRUE, stringsAsFactors=FALSE))
fnames <- gsub("topTable_|.txt", "", basename(top_files))
names(top) <- fnames

# Annotation file
annotation <- gmtPathways(annotation_file)

# Keep only GOBP + Reactome pathways
to_keep <- c("%GOBP|%Reactome")
annotation <- annotation[grep(to_keep, names(annotation), ignore.case=TRUE)]

# Filter out pathways with <10 genes and >500
annotation_length <- lapply(annotation, length)
annotation <- annotation[-which(annotation_length < path_min | annotation_length > path_max)]

######
# WORK BEGINS
######

# Sink log file to output folder
sink_file <- sprintf("%s/runEnrich.log", output_folder)
sink(sink_file)

# Date/time info
cat("RUNNING GENE ENRICHMENT/DEPLETION ANALYSIS\n")
Sys.time()

# Write enrichment parameters to log
cat("\nEnrichment analysis parameters:\n")
cat(sprintf(
 "\tInput directory: %s
  \tOutput directory: %s
  \tAnnotation file: %s
  \tPathway size limit: %s - %s genes
  \tLog2foldchange threshold: %s
  \tFDR threshold: %s
  \tAvg expression threshold: %s\n",
  main_folder, normalizePath(output_folder), annotation_file, path_min, path_max,
  lfc_cut, fdr_cut, exp_cut
))

#########################
#### ENRICHMENT ANALYSIS
#########################

# Set up enrichment function
cat("\n\nStarting up...\n\n")
runEnrich <- function(set_name) {

  # Define gene expression set
  data <- top[[set_name]]

  # Get background set of genes
  background <- data$EntrezGene.ID

  # Split into up/down regulated gene sets
  genes_up <- data %>%
    filter(logFC > lfc_cut, AveExpr > 0, adj.P.Val < fdr_cut) %>%
    select(EntrezGene.ID) %>%
    unlist %>%
    as.character

  genes_down <- data %>%
    filter(logFC < -lfc_cut, AveExpr > 0, adj.P.Val < fdr_cut) %>%
    select(EntrezGene.ID) %>%
    unlist %>%
    as.character

  # Put genes into list
  genes <- list()
  genes[["upregulated"]] <- genes_up
  genes[["downregulated"]] <- genes_down

  res_list <- list()
  for (k in seq_along(genes)) {

    # Define test genes in set
    set_genes <- genes[[k]]

    # Define output data sheet properties
    sheet_name <- names(genes[k])
    cat(sprintf("Expression set run: %s, %s\n", set_name, sheet_name))

    if (length(set_genes) == 0) {
      cat(" ** no genes remaining to calculate enrichment/depletion **\n")
      res_list[[sheet_name]] <- "N/A"
      next
    }

    # Prepare results table
    res <- data.frame(
      pathway=NA, size=NA, expected_total=NA, expected_num=NA, expected_frac=NA,
      real_total=NA, real_num=NA, real_frac=NA, enrichment=NA, pval=NA, FDR=NA, genes=NA
    )

    # Run for each pathway
    for (i in seq_along(annotation)) {

      # Define pathway being tested
      annotation_i <- annotation[i]
      annotation_name <- names(annotation_i)
      annotation_gene <- annotation[[i]]

      # Number of test genes in pathway
      real <- intersect(set_genes, annotation_gene) # overlap of genes in pathway
      real_frac <- (length(real) / length(set_genes)) * 100

      # Number of universe genes in pathway
      exp <- intersect(background, annotation_gene) # overlap of universe in pathway
      exp_frac <- (length(exp) / length(background)) * 100

      # Indicate whether term is enriched/depleted for genes
      if (real_frac > exp_frac) {
        enrich <- "enriched"
      } else if (exp_frac > real_frac) {
        enrich <- "depleted"
      } else if (exp_frac == real_frac) {
        enrich <- "none"
      } else {
        cat("what happened here ... ")
      }

      # Set up contingency table for p val calculation
      mat <- matrix(ncol=2, nrow=2)
      mat[1,] <- c(length(real), length(exp))
      mat[2,] <- c(length(set_genes), length(background))

      # P val
      pval <- fisher.test(mat, alternative="two.sided")

      # Vector of test genes in pathway
      if (length(real) == 0) {
        pathway_genes <- NA
      } else if (length(real) == 1) {
        pathway_genes <- real
      } else {
        pathway_genes <- paste(real, collapse="; ")
      }

      # Fill data in table
      res[i,]$expected_total <- length(background)
      res[i,]$real_total <- length(set_genes)
      res[i,]$pathway <- annotation_name
      res[i,]$size <- length(annotation_gene)
      res[i,]$expected_num <- length(exp)
      res[i,]$expected_frac <- round(exp_frac, 2)
      res[i,]$real_num <- length(real)
      res[i,]$real_frac <- round(real_frac, 2)
      res[i,]$enrichment <- enrich
      res[i,]$pval <- pval$p.value # do not round!
      res[i,]$genes <- pathway_genes
    }

    # FDR correct p values
    res <- res[order(res$pval),]
    res$FDR <- p.adjust(res$pval, method="BH") # do not round!

    # Check if pathways of interest come up, and print message to screen if so
    pathway_interest <- c("Irf|Nfk|Nfe|Nfr")
    pathways_found <- res[grep(pathway_interest, res$pathway, ignore.case=TRUE),]

    # Are the pathways enriched or depleted?
    n_sig_enrich <- length(which(pathways_found$enrichment == "enriched" & pathways_found$FDR < fdr_cut))
    if (n_sig_enrich != 0) {
      cat(sprintf("** %s Irf/Nfk/Nfe/Nfr pathway(s) significantly enriched **\n", n_sig_enrich))
    }

    n_sig_deplete <- length(which(pathways_found$enrichment == "depleted" & pathways_found$FDR < fdr_cut))
    if (n_sig_deplete != 0) {
      cat(sprintf("** %s Irf/Nfk/Nfe/Nfr pathway(s) significantly depleted **\n", n_sig_deplete))
    }

    # Write out to list
    res_list[[set_name]][[sheet_name]] <- res
  }
  return(res_list)
}

# Run enrichment on all combos of data
res_list_all <- mapply(runEnrich, fnames[17])

# Write out
cat("\n\nWriting excel files with results...\n")
res_fname <- sprintf("%s/enrich_%s_mouse_GOBP_reactome.xlsx", output_folder, fnames)
for (i in seq_along(res_fname)) {
  print(res_fname[i])
  write.xlsx(res_list_all[[1]], file=res_fname[i])
}

# End sinking
sink()
