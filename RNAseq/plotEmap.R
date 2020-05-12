#!/usr/bin/env Rscript
# Create EnrichmentMap in Cytoscape
# NOTE need to run runManualEnrich.R and writeEmapFile.R beforehand

# Loads all packages in a way that allows exporting to child environments
packages <- c("RCy3", "openxlsx")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only  =  TRUE))
}

source("../../usefulFunctions.R")

# Set working directory
setwd("/Users/catherineross/Dropbox/moffat/KL2020_analysis_CR/RNAseq")

######
# DEFINE I/O
######

annotation_file <- "input/pathway_annotations/Mouse_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt"
em_files <- list.files(pattern="table", path="output/em_files", full.names=TRUE)

######
# PREPARE DATA
######

# Set path to EM tables (enrichment results formatted for Cytoscape input)
em_files <- normalizePath(em_files)
em_names <- gsub("table_|.txt", "", basename(em_files))

######
# WORK BEGINS
######

# Confirm that Cytoscape is installed and opened
cytoscapePing()

# Initiate df of all pathway groupings
EM_group_df_all <- NULL
for (file in seq_along(em_files)) {

	# Map for each sample analysis
	EM_name <- em_names[file]
	EM_file <- em_files[file]
	print(EM_name)

	# Create EM using given parameters
	cat("* Building the network\n")
	em_command <- paste('enrichmentmap build analysisType="generic"',
			'gmtFile=', annotation_file,
			'enrichmentsDataset1=', EM_file,
			'pvalue=', 1,
			'qvalue=', 1,
			'similaritycutoff=', 0.375,
			'coefficients=', "COMBINED",
		  'combinedConstant=', 0.5)
	response <- commandsGET(em_command)

	# Annotate the network using AutoAnnotate app
	cat("* Annotating the network using AutoAnnotate\n")
	aa_command <- paste("autoannotate annotate-clusterBoosted",
		"clusterAlgorithm=MCL",
		"maxWords=3")
	#print(aa_command)
	response <- commandsGET(aa_command)

	# Pull node names from AutoAnnotated clusters
	# Get node and edge table that have the cluster annotations
	node_table <- getTableColumns(table="node", network=getNetworkSuid())
	edge_table <- getTableColumns(table="edge", network=getNetworkSuid())

	# Get the correct attribute names
	descr_attrib <- colnames(node_table)[grep(colnames(node_table), pattern="GS_DESCR")]

	# Get all the cluster numbers
	cluster_num <- node_table$'__mclCluster'

	# Get the unique set of clusters
	set_clusters <- unique(cluster_num)
	set_clusters <- set_clusters[which(set_clusters != 0)]

	# Calculate the cluster names using AutoAnnotate and collate all the nodes
	# associated with each cluster
	cluster_names <- c()
	for (i in 1:length(set_clusters)) {
	  current_cluster <- set_clusters[i]
	  gs_cluster <- node_table$name[which(node_table$'__mclCluster' == current_cluster)]

	  # For this cluster of gs get the gs descr to use in defining in AutoAnnotate
	  gs_cluster_suid <- node_table$SUID[which(node_table$name %in% gs_cluster)]
	  suids_aa <- paste("SUID", gs_cluster_suid, sep=":")

	  # Get the annotation for the given cluster
	  curernt_name = NULL
	  aa_label <- paste("autoannotate label-clusterBoosted",
				"labelColumn=", descr_attrib,
				"maxWords=3",
				"nodeList=", paste(suids_aa, collapse=","))

	  current_name <- commandsGET(aa_label)
	  cluster_names <- rbind(cluster_names,
			c(current_cluster, current_name, length(gs_cluster), paste(gs_cluster, collapse=";"))
		)
	}

	if (unique(cluster_names[,2]) != "") {
		# Create list of cluster names and associated pathways
		EM_group_list <- list()
		for (k in 1:nrow(cluster_names)) {
			EM_group_list[[k]] <- unlist(strsplit(cluster_names[k,4], ";"))
			names(EM_group_list)[k] <- gsub(" ", "_", toupper(cluster_names[k,2]))
		}
		# Get single unclustered gene sets and add to list
		single_path <- node_table[-which(node_table$name %in% unlist(EM_group_list)),]$name
	} else {
		EM_group_list <- NULL
		single_path <- node_table$name
	}

	# Remove annotation info from gene set names
	single_path_name <- gsub("\\%.*", "", single_path)
	single_path_name <- gsub(" ", "_", single_path_name)
	names(single_path) <- single_path_name

	# Add to existing list
	EM_group_list <- c(EM_group_list, single_path)
	EM_group_df <- as.data.frame(names(EM_group_list))
	colnames(EM_group_df) <- EM_name

	# Write out to master df
	EM_group_df_all <- cbind.na(EM_group_df_all, EM_group_df)
}

# Prepare for output
EM_group_df_all <- EM_group_df_all[,-1]
output_file <- "output/table_enrichment_EMgroups.xlsx"
write.xlsx(EM_group_df_all, output_file, row.names=FALSE)
