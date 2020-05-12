#!/usr/bin/env Rscript
# Analysis of in vivo mouse validation data

# Loads all packages in a way that allows exporting to child environments
packages <- c("openxlsx", "reshape2", "tidyr", "dplyr", "ggplot2", "ggforce", "cowplot", "RColorBrewer")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR/input")

######
# DEFINE I/O
######

## INPUTS ##

# In vivo data
info_file   <- "sample_info/samples_list_KL.xlsx"
rc_file     <- "readcounts/MUS018 EMT6-HA In Vivo Validation - 20200320 - readcounts.txt"
rc_val_file <- "readcounts/MUS014 EMT6-HA Validation - 20190409 - readcounts.txt"

# In vitro data
class_file   <- "library_annotations/Core_MeanNormZ.txt"
control_file <- "library_annotations/mVal_gene_list.txt"
bf_file      <- "library_annotations/bftable_all.txt"

## OUTPUTS ##

output_folder <- "/Users/catherineross/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVivo/unfiltered_data"
plot_folder   <- sprintf("%s/plots", output_folder)
table_folder  <- sprintf("%s/tables", output_folder)
fig_folder    <- sprintf("%s/clean_figs", output_folder)

######
# PARAMETER SETTING
######

# Sets special hard-coded column names
special_cols <- c("CHROMOSOME", "START", "STOP", "STRAND", "SEQUENCE", "GENE")

# Experimental control genes
control_genes <- c("luciferase", "LacZ", "EGFP")

# Timepoint groups
timepoints <- c("Early", "Late")

######
# DEFINE FUNCTIONS
######

# Log2 depth-normalize readcounts (KB)
normalizeReads <- function(data) {
  log2((sweep(data, 2, apply(data, 2, sum)/1000000, FUN = "/")) + 0.1)
}

# Detect outliers via PCA
# The standard way to detect outliers in genetics is the criterion of being
# “more than 6 standard deviations away from the mean”
# https://privefl.github.io/blog/detecting-outlier-samples-in-pca/
pcaOutlier <- function(df) {
  p_out <- apply(df, 2, function(x) which((abs(x - median(x)) / mad(x)) > 6)) %>% Reduce(union, .)
  return(p_out)
}

# Detect outliers in boxplots
boxOutlier <- function(x) {
  b_out <- x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
  return(b_out)
}

######
# READ DATA
######

# Library sample info
info <- read.xlsx(info_file)

# NOTE aggregating mid+late timepoint groups (now labelled "Late")
info$Group <- gsub("Mid", "Late", info$Group)

# In vivo readcounts
rc <- read.delim(rc_file, h = TRUE, stringsAsFactors = FALSE)

# In vitro readcounts
rc_val  <- read.delim(rc_val_file, h = TRUE, stringsAsFactors = FALSE)

# Subset in vivo library for those in validation library
rc <- rc[which(rc$SEQUENCE %in% rc_val$SEQUENCE),]

# In vitro bayes factor data
bf <- read.delim(bf_file, h = TRUE, stringsAsFactors = FALSE)

# Library annotation files
classes <- read.delim(class_file, h = TRUE, stringsAsFactors = FALSE)
controls <- read.delim(control_file, h = TRUE, stringsAsFactors = FALSE)

######
# PREAPRE MVAL LIBRARY ANNOTATIONS
######

# Complete validation library
gene_lib <- data.frame(GENE = unique(rc$GENE))

# Define suppressors and sensitizers using EMT6 parental line
## > 0 = supp; < 0 = sens
classes2 <- classes[,c("GENE", "EMT6")]
classes2 <- na.omit(classes2)
classes2$CLASS <- NA
classes2[which(classes2$EMT6 > 0), "CLASS"] <- "Suppressor"
classes2[which(classes2$EMT6 < 0), "CLASS"] <- "Sensitizer"
classes3 <- classes2[,-2]

# Controls category file
controls2 <- controls[,-c(1,4)] # remove core/intergenic classes
colnames(controls2) <- c("Targeting_control", "Non_targeting_control", "Other_control")

# Melt to tidy frame
controls_melt <- melt(controls2, measure.vars = colnames(controls2))

# Remove empty values
controls_melt[which(controls_melt$value == ""),] <- NA
controls_melt <- na.omit(controls_melt)
colnames(controls_melt) <- c("CLASS", "GENE")

# Fix gene names for indexing
controls_melt$GENE[controls_melt$GENE == "Luciferase"] <- "luciferase"
controls_melt$GENE[controls_melt$GENE == "CD274"] <- "Cd274"
controls_melt$GENE[controls_melt$GENE == "CD47"] <- "Cd47"
controls_melt <- controls_melt[,c("GENE", "CLASS")]

# Combine sens/supp + control dfs together
all_classes <- rbind(classes3, controls_melt)

# Merge with gene_lib df
gene_lib <- left_join(gene_lib, all_classes, by="GENE")

# Define genes without scores in EMT6 parental line as "Non_hit"
gene_lib[which(is.na(gene_lib$CLASS)), "CLASS"] <- "Non_hit"

######
# DEFINE ESSENTIALS
######

# Essential genes -- bayes factor > 50 in EMT6 parental line
bf_em <- bf[,c("GENE", "EMT6.HA_T11")]
bf_em_ess <- bf_em[which(bf_em[, "EMT6.HA_T11"] > 50), "GENE"]

# Add essentiality info to gene_lib df
gene_lib[which(gene_lib$GENE %in% bf_em_ess), "CLASS"] <- paste(gene_lib[which(gene_lib$GENE %in% bf_em_ess), "CLASS"], "(essential)")

# Define factor levels
class_levels <- c("Suppressor (essential)", "Suppressor",
                  "Sensitizer (essential)", "Sensitizer",
                  "Non_hit (essential)", "Non_hit",
                  "Targeting_control", "Non_targeting_control", "Other_control")

class_cols <- c("#F6EB13", "#999206", "#6D90CA", "#425A80",
                "#e62525", "#751d1d", "#C0BFBF", "forestgreen","black")

# Write out table
#table_out <- sprintf("%s/table_gene_classes_EMT6defined.xlsx", table_folder)
#write.xlsx(gene_lib, table_out, row.names = FALSE)

######
# SAMPLE LABELS
######

# Generate sample labels
groups <- info$Group
groups[grep("T", groups)] <- paste0(groups[grep("T", groups)], "_control")
replicates <- info$Replicate
conditions <- gsub(" ", "-", info$Condition)

# Put info together
labels <- paste(conditions, replicates, sep = "-")
labels <- paste(labels, groups, sep = "_")
labels[1:4] <- c("T0_control", "T6-A_control", "T6-B_control", "T6-C_control")
labels_df <- data.frame(sample = info$Sample, label = labels)
labels_df$label <- as.character(labels_df$label)

# Rename readcount sample labels
for (col in (length(special_cols)+1):ncol(rc)) {
  sample <- labels_df[which(labels_df$sample %in% colnames(rc)[col]), "label"]
  colnames(rc)[col] <- sample
}

# Merge T6 columns
rc[,"T6_control"] <- rowMeans(rc[,grep("T6", colnames(rc))], na.rm = TRUE)

######
# T6-NORMALIZED T0
######

# Define sample columns
t0_col <- "T0_control"
t6_col <- "T6_control"

# Recode missing values to 0
rc[is.na(rc)] <- 0

# Normalize T0 to T6 inputs
rc[,"T0-T6_control"] <- rc[[t0_col]]-rc[[t6_col]]

# Clean and re-order columns
rc2 <- rc[,-grep("T6-|vitro", colnames(rc))]
rc2 <- rc2[,c(seq_along(special_cols),7,186,187,8:185)]

# Define more sample columns
t0_t6_col <- "T0-T6_control"
control_cols <- c(t0_col, t6_col, "T0-T6_control")
balb_cols <- grep("BALB", colnames(rc2), value = TRUE)
nsg_cols <- grep("NSG", colnames(rc2), value = TRUE)

######
# RAW READCOUNT QC
######

# Checks readcount depth for all samples
cat(paste("T0 total reads:", sum(rc2[[t0_col]]), "\n"))
cat(paste("T6 total reads:", sum(rc2[[t6_col]]), "\n"))

for (col in nsg_cols) {
  total_reads <- sum(rc2[[col]])
  cat(paste(col, "total reads:", total_reads, "\n"))
}
for (col in balb_cols) {
  total_reads <- sum(rc2[[col]])
  cat(paste(col, "total reads:", total_reads, "\n"))
}

# Check representation across chromosomes
# Get number of reads per chromosome divided by number of guides
chrom_stats <- rc2[!(rc2$GENE %in% control_genes),] %>%
  group_by(CHROMOSOME) %>%
  summarize(total_reads = sum(!!as.name(col)), n = n())

chrom_stats <- chrom_stats[order(chrom_stats$n, decreasing = TRUE),]
chrom_stats$adj_reads <- chrom_stats$total_reads / chrom_stats$n

# Set factor levels
chrom_stats$CHROMOSOME <- factor(chrom_stats$CHROMOSOME, levels = rev(chrom_stats$CHROMOSOME))

# Generate table of guides per gene
guides_total <- as.data.frame(table(rc2$GENE))
colnames(guides_total) <- c("GENE", "N_total_guidess")

# Write out
#table_out <- sprintf("%s/table_gene_gRNA_representation_allSamples.xlsx", table_folder)
#write.xlsx(guides_total, table_out, row.names = FALSE)

######
# DATA PROCESSING
######

# Filter guides by t0 abundance
lower_threshold <- 40
upper_threshold <- 1e5
to_remove <- rep(FALSE, nrow(rc2))
to_remove[rc2[[t0_col]] < lower_threshold | rc2[[t0_col]] > upper_threshold] <- TRUE

# Normalizes readcounts to sequencing depth
rc2_norm <- rc2
rc2_norm[!(colnames(rc2_norm) %in% special_cols)] <- normalizeReads(rc2_norm[!(colnames(rc2_norm) %in% special_cols)])

# Remove guides that did not pass T0 readcount thresholds.
# Must be performed after read-depth normalization
removed_guides_ind <- which(to_remove)
rc2_norm <- rc2_norm[!to_remove,]
cat(paste("Excluded a total of", sum(to_remove), "guides for t0 representation\n"))

######
# PCA
######

# Prepare normalized reads for computing
rc2_t_norm <- t(rc2_norm[!(colnames(rc2_norm) %in% special_cols)])

# PCA analysis with prcomp
rc2_pca_norm <- prcomp(rc2_t_norm)

# Prepare for plotting
rc2_pca_norm_df <- as.data.frame(rc2_pca_norm$x)
rc2_pca_norm_df$label <- rownames(rc2_pca_norm_df)
rownames(rc2_pca_norm_df) <- NULL

# Further separate out sample labels
rc2_pca_norm_df2 <- separate(rc2_pca_norm_df, col = label, into = c("sample", "group"), sep = "_")

# Group into NSG / BALB sets
rc2_pca_norm_nsg <- rc2_pca_norm_df2[grep("NSG", rc2_pca_norm_df2$sample),]
rc2_pca_norm_nsg$group <- factor(rc2_pca_norm_nsg$group, levels = timepoints)
rc2_pca_norm_nsg$outlier <- "No"
rc2_pca_norm_balb <- rc2_pca_norm_df2[grep("BALB", rc2_pca_norm_df2$sample),]
rc2_pca_norm_balb$group <- factor(rc2_pca_norm_balb$group, levels = timepoints)
rc2_pca_norm_balb$outlier <- "No"

# Measure outlier-ness among each group of mice (eg NSG early, BALBc early, etc)
for (set in timepoints) {
  # NSG set
  nsg_set <- which(rc2_pca_norm_nsg$group == set)
  nsg_outliers <- pcaOutlier(rc2_pca_norm_nsg[nsg_set,1:2])
  rc2_pca_norm_nsg[nsg_set,][nsg_outliers, "outlier"] <- "Yes"
  # BALBc set
  balb_set <- which(rc2_pca_norm_balb$group == set)
  balb_outliers <- pcaOutlier(rc2_pca_norm_balb[balb_set,1:2])
  rc2_pca_norm_balb[balb_set,][balb_outliers, "outlier"] <- "Yes"
}

# List of sample outliers
nsg_outliers <- rc2_pca_norm_nsg[which(rc2_pca_norm_nsg$outlier == "Yes"), "sample"]
balb_outliers <- rc2_pca_norm_balb[which(rc2_pca_norm_balb$outlier == "Yes"), "sample"]
sample_outliers <- c(nsg_outliers, balb_outliers)

# Write out
#table_out <- sprintf("%s/list_pca_normReads_sample_outliers_NSG_BALB.txt", table_folder)
#write.table(sample_outliers, table_out)

######
# FILTER OUT "BAD" MICE
######

# Filter out only NSG mice (PCA defined)
#nsg_outliers <- grep("NSG", sample_outliers, value = TRUE)
#rc2_nsg_filter <- rc2[-grep(paste(nsg_outliers, collapse = "|"), colnames(rc2))]

# Filter out both NSG+BALBC mice (PCA defined)
#rc2_nsg_balb_filter <- rc2[-grep(paste(sample_outliers, collapse = "|"), colnames(rc2))]

######
# LOG2 NORMALIZED READS
######

## NSG ##
rc2_norm_nsg <- rc2_norm[which(colnames(rc2_norm) %in% c("GENE", "SEQUENCE", control_cols, nsg_cols))]
rc2_norm_nsg <- melt(rc2_norm_nsg)
rc2_norm_nsg <- separate(rc2_norm_nsg, col = variable, into = c("sample", "group"), sep = "_")
rc2_norm_nsg$group <- factor(rc2_norm_nsg$group, levels = c("control", timepoints))
rc2_norm_nsg <- left_join(rc2_norm_nsg, gene_lib)
rc2_norm_nsg$CLASS <- factor(rc2_norm_nsg$CLASS, levels = class_levels)

# Sample level gene mean
rc2_norm_nsg_mean <-
  rc2_norm_nsg %>%
  group_by(sample) %>%
  summarise(mean = median(value)) %>%
  as.data.frame()

# Order
rc2_norm_nsg_mean <- rc2_norm_nsg_mean[order(rc2_norm_nsg_mean$mean),]
rc2_norm_nsg$sample <- factor(rc2_norm_nsg$sample, levels = rc2_norm_nsg_mean$sample)

## BALB ##
rc2_norm_balb <- rc2_norm[which(colnames(rc2_norm) %in% c("GENE", "SEQUENCE", control_cols, balb_cols))]
rc2_norm_balb <- melt(rc2_norm_balb)
rc2_norm_balb <- separate(rc2_norm_balb, col = variable, into = c("sample", "group"), sep = "_")
rc2_norm_balb$group <- factor(rc2_norm_balb$group, levels = c("control", timepoints))
rc2_norm_balb <- left_join(rc2_norm_balb, gene_lib)
rc2_norm_balb$CLASS <- factor(rc2_norm_balb$CLASS, levels = class_levels)

# Sample level gene mean
rc2_norm_balb_mean <-
  rc2_norm_balb %>%
  group_by(sample) %>%
  summarise(mean = median(value)) %>%
  as.data.frame()

# Order
rc2_norm_balb_mean <- rc2_norm_balb_mean[order(rc2_norm_balb_mean$mean),]
rc2_norm_balb$sample <- factor(rc2_norm_balb$sample, levels = rc2_norm_balb_mean$sample)

######
# LOG2 FOLDCHANGE
######

# Calculate log2 foldchange
rc2_norm_t0 <- rc2_norm[[t0_col]]
rc2_lfc <- rc2_norm
rc2_lfc[!(colnames(rc2_lfc) %in% special_cols)] <- rc2_lfc[!(colnames(rc2_lfc) %in% special_cols)] - rc2_norm_t0

## NSG ##

# Subset for NSG/control columns
rc2_lfc_nsg <- rc2_lfc[which(colnames(rc2_lfc) %in% c("GENE", "SEQUENCE", nsg_cols))]

# Separate 'variable' column to extract timepoint groups
rc2_lfc_nsg2 <- melt(rc2_lfc_nsg)
rc2_lfc_nsg2 <- separate(rc2_lfc_nsg2, col = variable, into = c("sample", "group"), sep = "_")
rc2_lfc_nsg2$group <- factor(rc2_lfc_nsg2$group, levels = c("control", timepoints))
rc2_lfc_nsg2 <- left_join(rc2_lfc_nsg2, gene_lib)
rc2_lfc_nsg2$CLASS <- factor(rc2_lfc_nsg2$CLASS, levels = class_levels)

# Get means
rc2_lfc_nsg_mean <-
  rc2_lfc_nsg2 %>%
  group_by(GENE, group, sample) %>%
  summarise(mean = mean(value)) %>%
  as.data.frame()

rc2_lfc_nsg3 <- left_join(rc2_lfc_nsg2, rc2_lfc_nsg_mean)

## BALB C ##

# Subset for BALBC/control columns
rc2_lfc_balb <- rc2_lfc[which(colnames(rc2_lfc) %in% c("GENE", "SEQUENCE", balb_cols))]

# Separate 'variable' column to extract timepoint groups
rc2_lfc_balb2 <- melt(rc2_lfc_balb)
rc2_lfc_balb2 <- separate(rc2_lfc_balb2, col = variable, into = c("sample", "group"), sep = "_")
rc2_lfc_balb2$group <- factor(rc2_lfc_balb2$group, levels = c("control", timepoints))
rc2_lfc_balb2 <- left_join(rc2_lfc_balb2, gene_lib)
rc2_lfc_balb2$CLASS <- factor(rc2_lfc_balb2$CLASS, levels = class_levels)

# Get means
rc2_lfc_balb_mean <-
  rc2_lfc_balb2 %>%
  group_by(GENE, group, sample) %>%
  summarise(mean = mean(value)) %>%
  as.data.frame()

rc2_lfc_balb3 <- left_join(rc2_lfc_balb2, rc2_lfc_balb_mean)

######
# T6 NORMALIZED LOG2 FOLDCHANGE
######

# Calculate log2 foldchange using T6-normalized T0
rc2_norm_t0_t6 <- rc2_norm[[t0_t6_col]]
rc3_lfc <- rc2_norm
rc3_lfc[!(colnames(rc3_lfc) %in% special_cols)] <- rc3_lfc[!(colnames(rc3_lfc) %in% special_cols)] - rc2_norm_t0_t6

## NSG ##

# Subset for NSG/control columns
rc3_lfc_nsg <- rc3_lfc[which(colnames(rc3_lfc) %in% c("GENE", "SEQUENCE", nsg_cols))]

# Separate 'variable' column to extract timepoint groups
rc3_lfc_nsg2 <- melt(rc3_lfc_nsg)
rc3_lfc_nsg2 <- separate(rc3_lfc_nsg2, col = variable, into = c("sample", "group"), sep = "_")
rc3_lfc_nsg2$group <- factor(rc3_lfc_nsg2$group, levels = c("control", timepoints))
rc3_lfc_nsg2 <- left_join(rc3_lfc_nsg2, gene_lib)
rc3_lfc_nsg2$CLASS <- factor(rc3_lfc_nsg2$CLASS, levels = class_levels)

# Get means
rc3_lfc_nsg_mean <-
  rc3_lfc_nsg2 %>%
  group_by(GENE, group, sample) %>%
  summarise(mean = mean(value)) %>%
  as.data.frame()

rc3_lfc_nsg3 <- left_join(rc3_lfc_nsg2, rc3_lfc_nsg_mean)

## BALBC ##

# Subset for BALBC/control columns
rc3_lfc_balb <- rc3_lfc[which(colnames(rc3_lfc) %in% c("GENE", "SEQUENCE", balb_cols))]

# Separate 'variable' column to extract timepoint groups
rc3_lfc_balb2 <- melt(rc3_lfc_balb)
rc3_lfc_balb2 <- separate(rc3_lfc_balb2, col = variable, into = c("sample", "group"), sep = "_")
rc3_lfc_balb2$group <- factor(rc3_lfc_balb2$group, levels = c("control", timepoints))
rc3_lfc_balb2 <- left_join(rc3_lfc_balb2, gene_lib)
rc3_lfc_balb2$CLASS <- factor(rc3_lfc_balb2$CLASS, levels = class_levels)

# Get means
rc3_lfc_balb_mean <-
  rc3_lfc_balb2 %>%
  group_by(GENE, group, sample) %>%
  summarise(mean = mean(value)) %>%
  as.data.frame()

rc3_lfc_balb3 <- left_join(rc3_lfc_balb2, rc3_lfc_balb_mean)

######
# DIFFERENTIAL LFC (BALB - NSG)
######

## PER SAMPLE ##

# Construct df
rc3_lfc_diff <- data.frame(
  GENE = rc3_lfc_balb2$GENE,
  SEQUENCE = rc3_lfc_balb2$SEQUENCE,
  group = rc3_lfc_balb2$group,
  sample = rc3_lfc_balb2$sample,
  LFC_diff = rc3_lfc_balb2$value
)

# Mean LFC of NSG samples
rc3_lfc_nsg_mean <-
  rc3_lfc_nsg3 %>%
  group_by(GENE, group) %>%
  summarise(nsg_mean = mean(mean)) %>%
  as.data.frame()

# Mean LFC of NSG samples
rc3_lfc_balb_mean <-
  rc3_lfc_balb3 %>%
  group_by(GENE, group) %>%
  summarise(balb_mean = mean(mean)) %>%
  as.data.frame()

# Subtract BALBc LFC - mean NSG LFC
for (set in unique(rc3_lfc_diff$group)) {
  for (gene in unique(rc3_lfc_diff$GENE)) {
    nsg_mean_set <- rc3_lfc_nsg_mean[which(rc3_lfc_nsg_mean$group == set & rc3_lfc_nsg_mean$GENE == gene), "nsg_mean"]
    rc3_lfc_diff[which(rc3_lfc_diff$group == set & rc3_lfc_diff$GENE == gene), "LFC_diff"] <-
      (rc3_lfc_diff[which(rc3_lfc_diff$group == set & rc3_lfc_diff$GENE == gene), "LFC_diff"] - nsg_mean_set)
  }
}

# Get gene-level means per sample
gene_means <-
  rc3_lfc_diff %>%
  group_by(GENE, sample) %>%
  summarise(gene_mean = mean(LFC_diff)) %>%
  as.data.frame()

# Merge
rc3_lfc_diff2 <- left_join(rc3_lfc_diff, gene_means)

## PER TIMEPOINT GROUP ##

# SEQUENCE-level mean of sample diff LFC per timepoint group (aggregating samples)
rc3_lfc_diff_mean <-
  rc3_lfc_diff %>%
  group_by(GENE, SEQUENCE, group) %>%
  summarise(mean = mean(LFC_diff), sd = sd(LFC_diff)) %>%
  as.data.frame()

# Get gene-level means per group
gene_means2 <-
  rc3_lfc_diff_mean %>%
  group_by(GENE, group) %>%
  summarise(gene_mean = mean(mean)) %>%
  as.data.frame()

# Merge
rc3_lfc_diff_mean2 <- left_join(rc3_lfc_diff_mean, gene_means2)

# Get mean of non-targetting controls at each timepoint group
control_means2 <-
  rc3_lfc_diff_mean2 %>%
  filter(GENE == "luciferase" | GENE == "EGFP" | GENE == "LacZ") %>%
  group_by(group) %>%
  summarise(control_mean = mean(mean))

# Merge annotations with diffLFC dataframes
## SAMPLE LEVEL
rc3_lfc_diff3 <- left_join(rc3_lfc_diff2, gene_lib)
rc3_lfc_diff3$CLASS <- factor(rc3_lfc_diff3$CLASS, levels = class_levels)

## GROUP LEVEL
rc3_lfc_diff_mean3 <- left_join(rc3_lfc_diff_mean2, gene_lib)
rc3_lfc_diff_mean3$CLASS <- factor(rc3_lfc_diff_mean3$CLASS, levels = class_levels)

# Define outliers
rc3_lfc_diff_mean3 <-
  rc3_lfc_diff_mean3 %>%
  group_by(group, CLASS) %>%
  mutate(outlier = ifelse(boxOutlier(gene_mean), GENE, NA)) %>%
  as.data.frame()

rc3_lfc_diff_mean3[!is.na(rc3_lfc_diff_mean3$outlier), "outlier"] <- as.character(rc3_lfc_diff_mean3[!is.na(rc3_lfc_diff_mean3$outlier), "GENE"])

######
# PLOT
######

setwd("/Users/catherineross/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVivo/unfiltered_data/plots")

# 1) Distribution of NSG log2 normalized foldchange across samples
plot1_name <- "plot_boxplot_sample_log2FC_NSG_control.pdf"
pdf(plot1_name, width = 20, height = 5, useDingbats = FALSE)

ggplot(rc2_lfc_nsg3, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange", x = "Mouse sample",
      title = "Distribution of gRNA log2 foldchange across NSG mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 2) Distribution of NSG log2 + T6 normalized foldchange across samples
plot2_name <- "plot_boxplot_sample_log2FC_T6norm_NSG_control.pdf"
pdf(plot2_name, width = 15, height = 4, useDingbats = FALSE)

ggplot(rc3_lfc_nsg3, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange\n(normalized to T6)", x = "Mouse sample",
      title = "Distribution of T6-normalized gRNA log2 foldchange across NSG mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 3) Distribution of BALB-C log2 normalized foldchange across samples
plot3_name <- "plot_boxplot_sample_log2FC_BALBC_control.pdf"
pdf(plot3_name, width = 20, height = 5, useDingbats = FALSE)

ggplot(rc2_lfc_balb3, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange", x = "Mouse sample",
      title = "Distribution of gRNA log2 foldchange across BALB-C mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 4) Distribution of BALB-C log2 + T6 normalized foldchange across samples
plot4_name <- "plot_boxplot_sample_log2FC_T6norm_BALBC_control.pdf"
pdf(plot4_name, width = 15, height = 4, useDingbats = FALSE)

ggplot(rc3_lfc_balb3, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange\n(normalized to T6)", x = "Mouse sample",
      title = "Distribution of T6-normalized gRNA log2 foldchange across BALB-C mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 5) Long scatterplot of gene differential LFC (BALB-C - NSG)
plot5_name <- "plot_scatter_gene_difflog2FC_T6norm_BALBC_NSG_classCol.pdf"
pdf(plot5_name, width = 45, height = 5, useDingbats = FALSE)

# Plot separately per group
for (set in timepoints) {
  # Prepare data
  df_plot <- filter(rc3_lfc_diff_mean3, group == set)
  df_plot <- df_plot[order(df_plot$gene_mean),]
  df_plot$GENE <- factor(df_plot$GENE, levels = unique(as.character(df_plot$GENE)))
  # Plot
  p <- ggplot(df_plot, aes(x = GENE, y = mean, colour = CLASS)) +
        geom_point() +
        geom_hline(data = control_means2[which(control_means2$group==set),], aes(yintercept = control_mean), colour = "green", linetype = "dashed") +
        labs(y = "Differential log2 foldchange\n(normalized to T6)", x = "Library gene",
            title = sprintf("%s timepoint distribution of T6-normalized gRNA differential log2 foldchange (BALB-C vs. NSG) across gene library", set),
            subtitle = "Reference line: gene-level mean of non-targeting controls") +
        theme_bw() +
        scale_colour_manual(values = class_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
}

dev.off()

# 6) Long scatterplot of gene differential LFC (BALB-C - NSG) per sample
plot6_name <- "plot_scatter_gene_difflog2FC_T6norm_BALBC_NSG_classCol_perSample.pdf"
pdf(plot6_name, width = 45, height = 5, useDingbats = FALSE)

# Plot separately per group
for (k in unique(rc3_lfc_diff3$sample)) {
  print(k)
  # Prepare data
  df_plot <- filter(rc3_lfc_diff3, sample == k)
  df_plot <- df_plot[order(df_plot$gene_mean),]
  df_plot$GENE <- factor(df_plot$GENE, levels = unique(df_plot$GENE))
  # Plot
  p <- ggplot(df_plot, aes(x = GENE, y = LFC_diff, colour = CLASS)) +
        geom_point() +
        #geom_hline(data = control_means[which(control_means$group==set),], aes(yintercept = control_mean), colour = "green", linetype = "dashed") +
        labs(y = "Differential log2 foldchange\n(normalized to T6)", x = "Library gene",
            title = sprintf("%s (%s)", as.character(unique(df_plot$sample)), as.character(unique(df_plot$group))),
            subtitle = "Distribution of T6-normalized gRNA differential log2 foldchange (BALB-C vs. NSG) across gene library") +
        theme_bw() +
        scale_colour_manual(values = class_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
}

dev.off()

# 7) Scatterplot of T0 vs. T6 readcounts
plot7_name <- "plot_scatter_T0_vs_T6_reads.pdf"
pdf(plot7_name, width = 5, height = 5, useDingbats = FALSE)

plot(rc2$T0_control, rc2$T6_control)
plot(rc2_norm$T0_control, rc2_norm$T6_control)

dev.off()

# 9) PCA plots of NSG vs. BALB early/mid/late log2 normalized readcounts
plot9_name <- "plot_pca_depthNormreads_NSG_BALBC_timepoint_groups_outliers.pdf"
#plot9_name <- "plot_pca_log2reads_NSG_BALBC_timepoint_groups_outliers.pdf"
pdf(plot9_name, width = 15, height = 5)

## NSG
ggplot(rc2_pca_norm_nsg, aes(x = PC1, y = PC2, colour = outlier, label = sample)) +
  facet_wrap(.~group, scales = "free") +
  geom_point() +
  geom_text(size = 3, vjust = "inward", hjust = "inward") +
  scale_colour_manual(values = c("black", "red")) +
  ggtitle("PCA of log2-normalized readcounts across NSG mouse samples") +
  theme_bw() +
  theme(legend.position = "none")

## BALB
ggplot(rc2_pca_norm_balb, aes(x = PC1, y = PC2, colour = outlier, label = sample)) +
  facet_wrap(.~group, scales = "free") +
  geom_point() +
  geom_text(size = 3, vjust = "inward", hjust = "inward") +
  scale_colour_manual(values = c("black", "red")) +
  ggtitle("PCA of log2-normalized readcounts across BALB-C mouse samples") +
  theme_bw() +
  theme(legend.position = "none")

## BALB (without late 27-B and 10-A)
to_remove <- c("BALB/c-mouse-10-A", "BALB/c-mouse-27-B", "BALB/c-mouse-4-A",
               "BALB/c-mouse-33-B", "BALB/c-mouse-3-C")
rc2_pca_norm_balb2 <- rc2_pca_norm_balb[-grep(paste(to_remove, collapse = "|"), rc2_pca_norm_balb$sample),]

ggplot(rc2_pca_norm_balb2, aes(x = PC1, y = PC2, colour = outlier, label = sample)) +
  facet_wrap(.~group, scales = "free") +
  geom_point() +
  geom_text(size = 3, vjust = "inward", hjust = "inward") +
  scale_colour_manual(values = c("black", "red")) +
  labs(title = "PCA of depth-normalized readcounts across BALB-C mouse samples",
       subtitle = sprintf("Samples %s filtered out", paste(to_remove, collapse = ", "))) +
  theme_bw() +
  theme(legend.position = "none")

dev.off()

# 10) Boxplot of mean gene-level log2FC across timepoints
plot10_name <- "plot_boxplot_class_difflog2FC_T6norm_BALBC_NSG_classCol_outlier.pdf"
pdf(plot10_name, width = 12, height = 5)

ggplot(rc3_lfc_diff_mean3, aes(x = CLASS, y = gene_mean, fill = CLASS)) +
  facet_wrap(.~group, scales = "free") +
  geom_boxplot() +
  geom_text(aes(label = outlier), na.rm = TRUE, nudge_y = 0.1, nudge_x = 0.5) +
  geom_hline(data = control_means2, aes(yintercept = control_mean), colour = "green", linetype = "dashed") +
  labs(y = "Differential log2 foldchange\n(normalized to T6)",
      title = "Degree of suppressor and sensitizer T-cell dropout vs. controls",
      subtitle = "Reference line: gene-level mean of non-targeting controls") +
  theme_bw() +
  scale_fill_manual(values = class_cols) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

dev.off()

# 11) Boxplot of mean gene-level log2FC across timepoints per sample
plot11_name <- "plot_boxplot_class_difflog2FC_T6norm_BALBC_NSG_classCol_perSample.pdf"
pdf(plot11_name, width = 12, height = 10)

for (page in 1:8) {
  p <- ggplot(rc3_lfc_diff3, aes(x = CLASS, y = LFC_diff, fill = CLASS)) +
        facet_wrap_paginate(sample~group, scales = "free", ncol = 4, nrow = 3, page = page) +
        geom_boxplot() +
        labs(y = "Differential log2 foldchange\n(normalized to T6)",
            title = "Degree of suppressor and sensitizer T-cell dropout vs. controls") +
        theme_bw() +
        scale_fill_manual(values = class_cols) +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
}
dev.off()

# 12) Boxplot of log2 normalized readcounts for NSG mice
plot12_name <- "plot_boxplot_sample_log2reads_NSG_control.pdf"
pdf(plot12_name, width = 20, height = 5, useDingbats = FALSE)

ggplot(rc2_norm_nsg, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  labs(y = "Log2-normalized readcounts", x = "Mouse sample",
      title = "Distribution of log2-normalized reads across NSG mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 13) Boxplot of log2 normalized readcounts for BALBc mice
plot13_name <- "plot_boxplot_sample_log2reads_BALBC_control.pdf"
pdf(plot13_name, width = 20, height = 5, useDingbats = FALSE)

ggplot(rc2_norm_balb, aes(x = sample, y = value)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  labs(y = "Log2-normalized readcounts", x = "Mouse sample",
      title = "Distribution of log2-normalized reads across BALBc mouse samples at early/mid/late timepoints") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# 14) Rug-style scatterplot of log2 normalized readcounts per gene across NSG
plot14_name <- "plot_scatterRug_gene_log2reads_NSG_control.pdf"
pdf(plot14_name, height = 55, width = 5)

for (gene in unique(rc2_norm_nsg$GENE)) {
  df_plot <- filter(rc2_norm_nsg, GENE == gene)
  p <- ggplot(df_plot, aes(x = value, y = GENE)) +
        facet_wrap(.~sample, strip.position = "top", ncol = 1) +
        geom_point() +
        labs(x = "Log2-normalized readcounts",
             title = sprintf("%s (%s)", unique(df_plot$GENE), unique(df_plot$CLASS)),
             subtitle = "Distribution of guides across NSG mouse samples") +
        theme_bw()+
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
  print(p)
}
dev.off()

# 15) Rug-style scatterplot of log2 normalized readcounts per gene across BALB-C
plot15_name <- "plot_scatterRug_gene_log2reads_BALBC_control.pdf"
pdf(plot15_name, height = 55, width = 5)

for (gene in unique(rc2_norm_balb$GENE)) {
  df_plot <- filter(rc2_norm_balb, GENE == gene)
  p <- ggplot(df_plot, aes(x = value, y = GENE)) +
        facet_wrap(.~sample, strip.position = "top", ncol = 1) +
        geom_point() +
        labs(x = "Log2-normalized readcounts",
             title = sprintf("%s (%s)", unique(df_plot$GENE), unique(df_plot$CLASS)),
             subtitle = "Distribution of guides across BALB-C mouse samples") +
        theme_bw()+
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
  print(p)
}
dev.off()

# 16) Long scatterplot of log2 normalized reads per NSG sample
plot16_name <- "plot_scatter_gene_log2reads_NSG_classCol_perSample.pdf"
pdf(plot16_name, width = 45, height = 5, useDingbats = FALSE)

# Plot separately per group
for (k in unique(rc2_norm_nsg$sample)) {
  print(k)
  # Prepare data
  df_plot <- filter(rc2_norm_nsg, sample == k)
  df_plot <- df_plot[order(df_plot$value),]
  df_plot$GENE <- factor(df_plot$GENE, levels = unique(df_plot$GENE))
  # Plot
  p <- ggplot(df_plot, aes(x = GENE, y = value, colour = CLASS)) +
        geom_point() +
        #geom_hline(data = control_means[which(control_means$group==set),], aes(yintercept = control_mean), colour = "green", linetype = "dashed") +
        labs(y = "Log2-normalized reads", x = "Library gene",
            title = sprintf("%s (%s)", as.character(unique(df_plot$sample)), as.character(unique(df_plot$group))),
            subtitle = "Distribution of log2-normalized reads across gene library") +
        theme_bw() +
        scale_colour_manual(values = class_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
}

dev.off()

# 17) Long scatterplot of log2 normalized reads per NSG sample
plot17_name <- "plot_scatter_gene_log2reads_BALBC_classCol_perSample.pdf"
pdf(plot17_name, width = 45, height = 5, useDingbats = FALSE)

# Plot separately per group
for (k in unique(rc2_norm_balb$sample)) {
  print(k)
  # Prepare data
  df_plot <- filter(rc2_norm_balb, sample == k)
  df_plot <- df_plot[order(df_plot$value),]
  df_plot$GENE <- factor(df_plot$GENE, levels = unique(df_plot$GENE))
  # Plot
  p <- ggplot(df_plot, aes(x = GENE, y = value, colour = CLASS)) +
        geom_point() +
        #geom_hline(data = control_means[which(control_means$group==set),], aes(yintercept = control_mean), colour = "green", linetype = "dashed") +
        labs(y = "Log2-normalized reads", x = "Library gene",
            title = sprintf("%s (%s)", as.character(unique(df_plot$sample)), as.character(unique(df_plot$group))),
            subtitle = "Distribution of log2-normalized reads across gene library") +
        theme_bw() +
        scale_colour_manual(values = class_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.title = element_blank())
  print(p)
}

dev.off()

# 18) Boxplot distribution of NSG log2 + T6 normalized foldchange
plot18_name <- "plot_boxplot_class_log2FC_T6norm_NSG_classCol.pdf"
pdf(plot18_name, width = 10, height = 4, useDingbats = FALSE)

ggplot(rc3_lfc_nsg3, aes(x = CLASS, y = mean, fill = CLASS)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange\n(normalized to T6)",
       title = "Mean dropout of sensitizers, suppressors, essentials across NSG samples") +
  theme_bw() +
  scale_fill_manual(values = class_cols) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

dev.off()

# 19) Boxplot distribution of BALBC log2 + T6 normalized foldchange
plot19_name <- "plot_boxplot_class_log2FC_T6norm_BALBC_classCol.pdf"
pdf(plot19_name, width = 10, height = 4, useDingbats = FALSE)

ggplot(rc3_lfc_balb3, aes(x = CLASS, y = mean, fill = CLASS)) +
  facet_grid(.~group, scales = "free_x", space = "free_x") +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
  labs(y = "Log2 foldchange\n(normalized to T6)",
       title = "Mean dropout of sensitizers, suppressors, essentials across BALB-C samples") +
  theme_bw() +
  scale_fill_manual(values = class_cols) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

dev.off()

# 20)  Boxplot distribution of NSG log2 + T6 normalized foldchange per sample
plot20_name <- "plot_boxplot_class_log2FC_T6norm_NSG_classCol_perSample.pdf"
pdf(plot20_name, width = 15, height = 10, useDingbats = FALSE)

for (page in 1:8) {
  p <- ggplot(rc3_lfc_nsg3, aes(x = CLASS, y = mean, fill = CLASS)) +
          facet_wrap_paginate(sample~group, scales = "free", ncol = 4, nrow = 3, page = page) +
          geom_boxplot() +
          geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
          labs(y = "Log2 foldchange\n(normalized to T6)",
               title = "Mean dropout of sensitizers, suppressors, essentials across NSG samples") +
          theme_bw() +
          scale_fill_manual(values = class_cols) +
          theme(axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank())
  print(p)
}
dev.off()

# 21) Boxplot distribution of BALB-C log2 + T6 normalized foldchange per sample
plot21_name <- "plot_boxplot_class_log2FC_T6norm_BALBC_classCol_perSample.pdf"
pdf(plot21_name, width = 15, height = 10, useDingbats = FALSE)

for (page in 1:8) {
  p <- ggplot(rc3_lfc_balb3, aes(x = CLASS, y = mean, fill = CLASS)) +
          facet_wrap_paginate(sample~group, scales = "free", ncol = 4, nrow = 3, page = page) +
          geom_boxplot() +
          geom_hline(yintercept = 0, linetype = "dashed", colour = "green") +
          labs(y = "Log2 foldchange\n(normalized to T6)",
               title = "Mean dropout of sensitizers, suppressors, essentials across BALB-C samples") +
          theme_bw() +
          scale_fill_manual(values = class_cols) +
          theme(axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank())
  print(p)
}
dev.off()


# 22) Density rug plot of NSG log2 readcounts per sample
plot22_name <- "plot_densityRug_sample_log2reads_NSG_classCol_perSample_essentials.pdf"
pdf(plot22_name, width = 6, height = 100, useDingbats = FALSE)

essentials <- rc2_norm_nsg[grep("essential", rc2_norm_nsg$CLASS),]
ggplot(rc2_norm_nsg, aes(x = value)) +
  facet_wrap(.~sample, strip.position = "top", ncol = 1, scales = "free_y") +
  geom_density() +
  labs(x = "Log2-normalized readcounts", y = "Density",
       title = "Distribution of essential gene guides across NSG samples") +
  geom_rug(data = essentials, aes(x = value, group = CLASS, colour = CLASS), sides = "b", length = unit(0.1, "npc")) +
  scale_colour_manual(values = class_cols[c(1,3,5)]) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

dev.off()

# 23) Density rug plot of BALBC log2 readcounts per sample
plot23_name <- "plot_densityRug_sample_log2reads_BALBC_classCol_perSample_essentials.pdf"
pdf(plot23_name, width = 6, height = 100, useDingbats = FALSE)

essentials <- rc2_norm_balb[grep("essential", rc2_norm_balb$CLASS),]
ggplot(rc2_norm_balb, aes(x = value)) +
  facet_wrap(.~sample, strip.position = "top", ncol = 1, scales = "free_y") +
  geom_density() +
  labs(x = "Log2-normalized readcounts", y = "Density",
       title = "Distribution of essential gene guides across BALB-C samples") +
  geom_rug(data = essentials, aes(x = value, group = CLASS, colour = CLASS), sides = "b", length = unit(0.1, "npc")) +
  scale_colour_manual(values = class_cols[c(1,3,5)]) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

dev.off()


# 24) Density rug plot of NSG log2 foldchange per sample
plot24_name <- "plot_densityRug_sample_log2FC_T6norm_NSG_classCol_perSample_essentials.pdf"
pdf(plot24_name, width = 6, height = 100, useDingbats = FALSE)

essentials <- rc3_lfc_nsg3[grep("essential", rc3_lfc_nsg3$CLASS),]
ggplot(rc3_lfc_nsg3, aes(x = value)) +
  facet_wrap(.~sample, strip.position = "top", ncol = 1, scales = "free_y") +
  geom_density() +
  labs(x = "Log2 foldchange\n(normalized to T6)", y = "Density",
       title = "Distribution of essential gene guides across NSG samples") +
  geom_rug(data = essentials, aes(x = value, group = CLASS, colour = CLASS), sides = "b", length = unit(0.1, "npc")) +
  scale_colour_manual(values = class_cols[c(1,3,5)]) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

dev.off()

# 25) Density rug plot of NSG log2 foldchange per sample
plot25_name <- "plot_densityRug_sample_log2FC_T6norm_BALBC_classCol_perSample_essentials.pdf"
pdf(plot25_name, width = 6, height = 100, useDingbats = FALSE)

essentials <- rc3_lfc_balb3[grep("essential", rc3_lfc_balb3$CLASS),]
ggplot(rc3_lfc_balb3, aes(x = value)) +
  facet_wrap(.~sample, strip.position = "top", ncol = 1, scales = "free_y") +
  geom_density() +
  labs(x = "Log2 foldchange\n(normalized to T6)", y = "Density",
       title = "Distribution of essential gene guides across BALB-C samples") +
  geom_rug(data = essentials, aes(x = value, group = CLASS, colour = CLASS), sides = "b", length = unit(0.1, "npc")) +
  scale_colour_manual(values = class_cols[c(1,3,5)]) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

dev.off()

######
# CLEANED SUPPLEMENTAL TABLES
######

# Join mean BALB LFC + mean NSG LFC dfs
rc3_lfc_bn_mean <- left_join(rc3_lfc_balb_mean, rc3_lfc_nsg_mean)
rc3_lfc_bn_mean$deltaLFC <- rc3_lfc_bn_mean$balb_mean - rc3_lfc_bn_mean$nsg_mean
rc3_lfc_bn_mean$pval <- NULL
rc3_lfc_bn_mean$fdr <- NULL

# Get p values
pval_list <- list()
for (set in timepoints) {

  # Separate by time point group
  rc3_set <- filter(rc3_lfc_diff_mean3, group == set)
  rc3_set_control <- filter(rc3_set, CLASS == "Targeting_control")

  # Get p values
  pvals <- data.frame(GENE = unique(rc3_set$GENE), pval = NA, fdr = NA)

  for (i in pvals$GENE) {
    gene_i <- filter(rc3_set, GENE == i)
    pval <- wilcox.test(gene_i$mean, rc3_set_control$mean, alternative = "l")$p.value
    pvals[which(pvals$GENE == i), "pval"] <- pval
  }

  # FDR correct pvalues
  pvals <- pvals[order(pvals$pval),]
  pvals$fdr <- p.adjust(pvals$pval, method = "BH")
  pvals$fdr <- signif(pvals$fdr, 1)

  # Write out to list
  pval_list[[set]] <- pvals
}

# Separate data by timepoint group
early <- filter(rc3_lfc_bn_mean, group == "Early")
early <- early[,-2]
early <- left_join(early, pval_list[["Early"]])
colnames(early) <- c("Gene", "Mean LFC BALB-C (early)", "Mean LFC NSG (early)", "DeltaLFC (early)", "Pvalue (early)", "FDR (early)")
late <- filter(rc3_lfc_bn_mean, group == "Late")
late <- late[,-2]
late <- left_join(late, pval_list[["Late"]])
colnames(late) <- c("Gene", "Mean LFC BALB-C (late)", "Mean LFC NSG (late)", "DeltaLFC (late)", "Pvalue (late)", "FDR (late)")

# Combine together
rc3_lfc_bn_mean2 <- left_join(early, late)

# Write out
stable_out <- sprintf("%s/table_summary_gene_mean_LFC_BALBC_NSG_pwilcox_less_v3.xlsx", table_folder)
write.xlsx(rc3_lfc_bn_mean2, stable_out, row.names = FALSE)

######
# CLEANED PLOTS
######

setwd("/Users/catherineross/Dropbox/moffat/CTL_paper/Revised Manuscript/New figures/New data/Catherine/inVivo/unfiltered_data/clean_figs")

####
# Distribution of negative control deltaLFC (lacz, luciferase, egfp)
# SEQUENCE rug plots for all Atg genes at the bottom (separate early, mid, late)
####

# All 3 control sets in same density plot
control_dlfc_all <- filter(rc3_lfc_diff_mean3, CLASS == "Targeting_control")

# All genes of interest (for p val calculation)
all_genes <- unique(rc3_lfc_diff_mean3$GENE)

# Function to generate plots for 3 timepoints
plot_densityRug <- function(set, gene, width, height, fontsize) {

  # Define genes of interest to plot
  if (gene == "Atg") {
    #col_set <- "#F6EB13"
    gene_levels <- c("Atg3", "Atg7", "Atg5", "Atg10", "Atg9a", "Atg101", "Atg12", "Atg14",
                     "EGFP", "LacZ", "luciferase", "Controls")
    gene_levels2 <- gene_levels[-length(gene_levels)]
  }
  if (gene == "Fitm2") {
    #col_set <- "#F6EB13"
    gene_levels <- c("Fitm2", "EGFP", "LacZ", "luciferase", "Control")
    gene_levels2 <- gene_levels[-length(gene_levels)]
  }
  if (gene == "Socs1") {
    #col_set <- "#F6EB13"
    gene_levels <- c("Socs1", "EGFP", "LacZ", "luciferase", "Control")
    gene_levels2 <- gene_levels[-length(gene_levels)]
  }
  if (gene == "Suppressors") {
    #col_set <- "#6D90CA"
    gene_levels <- c("B2m", "Tap1", "Tap2", "Tapbp", "Jak1", "Jak2", "Stat1", "Stat2",
                     "Ifngr1", "Ifngr2", "EGFP", "LacZ", "luciferase", "Controls")
    gene_levels2 <- gene_levels[-length(gene_levels)]
  }
  if (gene == "other") {
    #col_set <- "#6D90CA"
    gene_levels <- c("Vdac2", "Ptpn2", "EGFP", "LacZ", "luciferase", "Control")
    gene_levels2 <- gene_levels[-length(gene_levels)]
  }

  # Plot data
  plot_data <- filter(rc3_lfc_diff_mean3, group == set)

  # Control data
  control_dlfc <- filter(plot_data, CLASS == "Targeting_control")
  control_dlfc$GENE <- "Controls"

  # Atg data
  gene_dlfc <- plot_data[grep(paste(gene_levels2, collapse = "|"), plot_data$GENE),]
  gene_dlfc$GENE <- factor(gene_dlfc$GENE, levels = gene_levels2)

  # Combine
  dlfc <- gene_dlfc

  for (i in gene_levels2) {
    control_dlfc_i <- control_dlfc
    control_dlfc_i[which(control_dlfc_i$GENE == "Controls"), "GENE"] <- i
    dlfc <- rbind(dlfc, control_dlfc_i)
  }

  # Significance values (diff of Atg vs control means)
  pvals <- data.frame(GENE = all_genes, pval = NA, fdr = NA)

  for (i in pvals$GENE) {
    gene_i <- filter(plot_data, GENE == i)
    pval <- wilcox.test(gene_i$mean, control_dlfc$mean, alternative = "l")$p.value
    pvals[which(pvals$GENE == i), "pval"] <- pval
  }

  # FDR correct pvalues
  pvals <- pvals[order(pvals$pval),]
  pvals$fdr <- p.adjust(pvals$pval, method = "BH")
  pvals$fdr <- signif(pvals$fdr, 1)

  # Combine with gene set of interest
  pvals <- pvals[which(pvals$GENE %in% dlfc$GENE),]

  # Order pval df by mid gene levels
  pvals <- pvals[order(match(pvals$GENE, gene_levels2)),]
  pvals$GENE <- factor(pvals$GENE, levels = gene_levels2)

  # Order genes
  dlfc$GENE <- factor(dlfc$GENE, levels = gene_levels2)

  # Group colours
  cols <- brewer.pal(8, "Set1")
  #cols <- cols[c(3,1,2)]
  cols <- cols[c(3,2)]

  # Median per group
  plot_means <- control_dlfc_all %>%
    group_by(group) %>%
    summarise(group_mean = mean(gene_mean))

  # X axis limits
  xlims = c(-14, 6)

  # Density plot all time point groups
  p1 <- ggplot(control_dlfc_all, aes(x = mean, colour = group, fill = group)) +
          geom_density(bw = 0.5, alpha = 0.5) +
          theme_classic(base_size = fontsize) +
          scale_colour_manual(values = cols) +
          scale_fill_manual(values = cols) +
          scale_x_continuous(limits = c(xlims[1], xlims[2]), breaks = seq(xlims[1], xlims[2], by = 2)) +
          labs(y = "Density", x = expression(paste(Delta, " log"[2], "(foldchange)"))) +
          geom_vline(data = plot_means[which(plot_means$group==set),], aes(xintercept = group_mean), size = 0.3, colour = "black", linetype = "dashed") +
          theme(legend.position = c(1.05, 0.7),
                legend.title = element_blank(),
                legend.text = element_text(size = fontsize/1.2),
                legend.key.size = unit(0.25, "line"),
                panel.grid = element_blank())

  if (set == "Early") { col_set <- cols[1] }
  if (set == "Late") { col_set <- cols[2] }

  # Rug plots
  p2 <- ggplot(gene_dlfc, aes(x = mean, colour = CLASS)) +
          facet_grid(GENE~., switch = "y") +
          theme_bw(base_size = fontsize) +
          scale_x_continuous(limits = c(xlims[1], xlims[2]), breaks = seq(xlims[1], xlims[2], by = 2)) +
          geom_rug(data = subset(dlfc, CLASS == "Targeting_control"), alpha = 0.1, col = "grey70", size = 0.7, length = unit(1, "npc")) +
          geom_rug(data = subset(dlfc, CLASS == "Non_targeting_control"), alpha = 1, col = col_set, size = 0.7, length = unit(1, "npc")) +
          geom_rug(data = subset(dlfc, CLASS != "Targeting_control" & CLASS != "Non_targeting_control"), alpha = 1, col = col_set, size = 0.7, length = unit(1, "npc")) +
          geom_vline(data = plot_means[which(plot_means$group==set),], aes(xintercept = group_mean), size = 0.3, colour = "black", linetype = "dashed") +
          theme(strip.text.y.left = element_text(angle = 0, face = "italic"),
                strip.background = element_rect(fill = "grey20", color = "grey20", size = 1),
                strip.text = element_text(colour = "white"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                panel.grid = element_blank(),
                legend.title = element_blank(),
                legend.position = "none")

  # Significance values
  gene_labs <- as.character(pvals$fdr)
  names(gene_labs) <- pvals$GENE

  p3 <- ggplot(pvals) +
          theme_minimal(base_size = fontsize) +
          facet_grid(GENE~., labeller = labeller(GENE = gene_labs)) +
          theme(strip.text.y.right = element_text(angle = 0, hjust = 0))

  # Draw out
  plot_name <- sprintf("plot_densityRug_posControlDist_guides%s_difflog2FC_all_%s_pwilcox_v7.pdf", gene, set)
  pdf(plot_name, width = width, height = height, useDingbats = FALSE)

  # Align density to rug plot
  left <- plot_grid(p1, p2, labels = NULL, ncol = 1, align = "v", axis = "bl", rel_heights = c(2,3))
  # Align rug to FDR plot
  right <- plot_grid(NULL, p3, ncol = 1, rel_heights = c(2,3))
  # Put together
  print(plot_grid(left, right, ncol = 2, rel_widths = c(7,0.8)))

  dev.off()

  # Clean up dlfc table to write out
  table_out <- gene_dlfc
  table_out <- table_out[,-8]
  table_out <- table_out[-grep("control", table_out$CLASS),]
  table_out <- table_out[order(table_out$GENE, table_out$mean),]

  # Write out
  table_name <- sprintf("table_guides%s_difflog2FC_all_%s_pwilcox_v7.xlsx", gene, set)
  write.xlsx(table_out, file = table_name, row.names = FALSE)
}

set_list <- timepoints
mapply(plot_densityRug, set_list, gene = "Atg", width = 2, height = 3, fontsize = 5)
mapply(plot_densityRug, set_list, gene = "Suppressors", width = 2, height = 3.5, fontsize = 5)
mapply(plot_densityRug, set_list, gene = "Socs1", width = 2, height = 2, fontsize = 5)
mapply(plot_densityRug, set_list, gene = "Fitm2", width = 2, height = 2, fontsize = 5)
mapply(plot_densityRug, set_list, gene = "other", width = 2, height = 2, fontsize = 5)


###
## Ranked gene dlog2FC scatterplot
###

for (set in timepoints) {

  print(set)
  test <- filter(rc3_lfc_diff_mean3, group == set)
  test <- test[order(test$gene_mean),]
  test <- test[,-c(2,4,5,8)]
  test <- unique(test)
  test$GENE <- factor(test$GENE, levels = test$GENE)

  # Mean of targeting controls
  control_mean <- mean(subset(test, CLASS == "Targeting_control")$gene_mean)

  # Label top genes on each side
  test$LABEL <- ""
  test <- test[order(test$gene_mean),]
  test[1:15, "LABEL"] <- as.character(test[1:15, "GENE"])
  test <- test[order(test$gene_mean, decreasing = TRUE),]
  test[1:15, "LABEL"] <- as.character(test[1:15, "GENE"])
  test$CLASS <- gsub(" \\(essential\\)", "", test$CLASS)
  test[(grep("control|hit", test$CLASS)), "CLASS"] <- "Control"

  # Grey-out insiginificant genes
  sd <- sd(test$gene_mean)
  test[which(test$gene_mean > control_mean-sd & test$gene_mean < control_mean+sd), "CLASS"] <- "NS"

  # Set factor levels
  test$CLASS <- factor(test$CLASS, levels = c("Suppressor", "Sensitizer", "Control", "NS"))

  p <- ggplot(test, aes(x = as.numeric(GENE), y = gene_mean, fill = CLASS, label = LABEL)) +
        geom_point(data = subset(test, CLASS == "NS"), aes(y = gene_mean, x = as.numeric(GENE)), size = 1, col = "grey") +
        geom_point(data = subset(test, CLASS == "Control"), aes(y = gene_mean, x = as.numeric(GENE)), size = 1, col = "grey") +
        geom_point(data = subset(test, CLASS == "Sensitizer"), aes(y = gene_mean, x = as.numeric(GENE)), size = 2, shape = 21, alpha = 0.7) +
        geom_point(data = subset(test, CLASS == "Suppressor"), aes(y = gene_mean, x = as.numeric(GENE)), size = 2, shape = 21, alpha = 0.7) +
        geom_hline(yintercept = control_mean, linetype = "dotted", col = "black") +
        scale_y_continuous(limits = c(-12, 4), breaks = seq(-12, 4, by = 2)) +
        theme_bw(base_size = 12) +
    #    scale_colour_manual(values = c("#C0BFBF", "#C0BFBF", "#6D90CA", "#F6EB13")) +
        scale_fill_manual(values = c("#C0BFBF", "#C0BFBF", "#6D90CA", "#F6EB13")) +
        labs(y = expression(paste(Delta, " log2-foldchange")), x = "Gene Rank") +
        theme(legend.position = "none",
              legend.title = element_blank(),
              panel.grid = element_blank())

    if (set == "Mid") {
      p <- p + geom_text_repel(data = subset(test, gene_mean < -2.5 & LABEL ! =  ""), nudge_x = 50, nudge_y = 1,
                               segment.color = "grey50", size = 3, direction = "both", hjust = 1) +
               geom_text_repel(data = subset(test, gene_mean > -2.5 & LABEL ! =  ""), nudge_x = -10, nudge_y = -1,
                              segment.color = "grey50", size = 3, direction = "both", hjust = -1)
    }
    if (set == "Early") {
      p <- p + geom_text_repel(data = subset(test, gene_mean < -1 & LABEL ! =  ""), nudge_x = 50, nudge_y = 0,
                               segment.color = "grey50", size = 3, direction = "both", hjust = 1) +
               geom_text_repel(data = subset(test, gene_mean > -1 & LABEL ! =  ""), nudge_x = -10, nudge_y = 0,
                              segment.color = "grey50", size = 3, direction = "both", hjust = -1)
    }
    if (set == "Late") {
      p <- p + geom_text_repel(data = subset(test, gene_mean < -6 & LABEL ! =  ""), nudge_x = 50, nudge_y = 0,
                               segment.color = "grey50", size = 3, direction = "both", hjust = 1) +
               geom_text_repel(data = subset(test, gene_mean > -6 & LABEL ! =  ""), nudge_x = -10, nudge_y = 0,
                               segment.color = "grey50", size = 3, direction = "both", hjust = -1)
    }

  # Draw out
  plot_name <- sprintf("plot_scatterRank_gene_difflog2FC_%s.pdf", set)
  pdf(plot_name, width = 4, height = 4)
  print(p)
  dev.off()
}
