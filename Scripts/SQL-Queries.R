# Script Description:
# Script with SQL-based to queries to extract the splice junctions of the genes
# harboring the splicing-associated variants of interest.
# The queries are gene-specific.

# List of required packages
required_packages <- c("optparse","dplyr", "readxl", "stringr", "sqldf", "tidyr", "knitr")

# Check if each package is installed and load it if available
#for (package in required_packages) {
#	if (!requireNamespace(package, quietly = TRUE)) {
#		install.packages(package, dependencies = TRUE)	# Install the package if not installed
#	}
#}


library(optparse)
# Option parsing
option_list <- list(
	make_option(c("-p", "--cohort_project"), type="character", help="Set to BeatAML or TCGA"),
	make_option(c("-c", "--cohort_metadata"), type="character", help="Path to cohort metadata csv file"),
	make_option(c("-m", "--mapping_option"), type="character", default=NULL, help="Specify the reads mapping type."),
	make_option(c("-i", "--splice_junction"), type="character", help="Path to annotated splice junction collection"),
	make_option(c("-o", "--output_directory"), type="character", help="Path to the output directory, creates a folder called ExtractedSJ-Annotated were writes the splice junctions")
)

opt_parser <- OptionParser(
	description =paste0("R Script with SQL-based to queries to extract the splice junctions of the genes harboring the splicing-associated variants of interest.\n",
											"It requires the Cohort RNASeq Metadata CSV file containing the samples' names and index columns: 'sample_id' and 'INDEX'.\n",
											"Define the cohort of interest to extract the splice junctions of the genes harboring potential SAVs of interest: 'BeatAML' or 'TCGA' \n",
											"Processing steps:\n",
											"1. SQL-based query to extract a particular gene splice junctions (manually defined based on SJ annotation)\n",
											"2. Query results processing (function: mainQueryProcess):\n",
											"\t 2.1 Separate the read counts array by sample and use the sample_id as the column name (function: LecturesArrayToSamples) \n",
											"\t 2.2 Dataframe formating and sample annotation (function: TransposeAnnotator) \n",
											"\t 2.3 Save the results as a table (function: writeQuery) \n"),

	usage = "Rscript SQL-queries.R -p <cohort_project> -c <cohort_metadata> [-m <mapping_option>]-i <splice_junction> -o <output_directory>",

	option_list = option_list)
opt <- parse_args(opt_parser)



# Check if all three arguments are defined
#if (is.null(opt$cohort_project) || is.null(opt$cohort_metadata) || is.null(opt$splice_junction) || is.null(opt$output_directory)) {
#	cat("Please provide values for all four arguments: -p <cohort_project>	-c <cohort_metadata> -i <splice_junction> -o <output_directory> .\n")
#	cat("For help, run: Rscript SQL-queries.R -h\n")
#	q(save = "no")
#}

# Set working environment:
## Load libraries
library(dplyr)
library(readr)
library(stringr)
library(sqldf)
library(tidyr)
library(knitr)

#sessionInfo()

# Get the Cohort Metadata, the SJ collection path and the directory path to save the extracted Splice Junctions
cohort_project <- opt$cohort_project
csv_file <- opt$cohort_metadata
SJ_collection <- opt$splice_junction
mapping_option <- opt$mapping_option # can be null if not defined
output_directory <- opt$output_directory # Creates a folder called "ExtractedSJ-Annotated

CohortMetadata <- read_csv(csv_file)
CohortMetadata <- as.data.frame(CohortMetadata) # Not a tibble
colnames(CohortMetadata)

# Check the splice junction collection file
tryCatch({
	if (file.exists(SJ_collection)) {
		cat("Located splice junction collection.\n")
	} else {
		cat("Splice junction collection not found.\n")
		stop("Splice junction collection not found.")
	}
}, error = function(e) {
	cat("An error occurred:", conditionMessage(e), "\n")
	stop("Error: Execution halted.")
})


# Check and add "/" to the output directory path if needed
tryCatch({
	if (dir.exists(output_directory)){
		if (!grepl("/$", output_directory)){output_directory <- paste0(output_directory, "/")}
		}else{
			stop("Output directory does not exist.\n")
		}
	}, error = function(e) {
	cat("An error occurred:", conditionMessage(e), "\n")
	stop("Error: Execution halted.")
})


if (!grepl("/$", output_directory)) {
	dir.exists(directory_path)
	output_directory <- paste0(output_directory, "/")
}

## Check and create "ExtractedSJ-Annotated" folder if it does not exist
annotatedSJ_folder <- file.path(output_directory, "ExtractedSJ-Annotated")
if (!file.exists(annotatedSJ_folder) || !file.info(annotatedSJ_folder)$isdir) {
	dir.create(annotatedSJ_folder, recursive = TRUE)
	cat("The 'ExtractedSJ-Annotated' folder has been created.\n")
} else {
	cat("The 'ExtractedSJ-Annotated' folder already exists.\n")
}
output_directory <- paste0(output_directory,"ExtractedSJ-Annotated/")

## Check the cohort_project:
if (cohort_project =="BeatAML"){
	cat("Extracting BeatAML SAV-related genes\n")
}else if(cohort_project == "TCGA"){
	cat("Extracting TCGA SAV-related genes\n")
}else{
	cat("No project selected\n")
}

cat(paste0("Defined arguments:\n\n"))
cat(paste0("Cohort project: ",cohort_project,"\n"))
cat(paste0("Splice junction collection: ",SJ_collection,"\n"))
cat(paste0("Output directory: ",output_directory,"\n\n"))


#### PROCESSING FUNCTIONS ####

## Selected annotations:
selectedAnnotations<- function(query_df){
	cat("\n")
	cat("Splice junction start gene annotations:\n")
	cat(unique(query_df$Gene_1))
	cat("\n")
	cat("Splice junction end gene annotations:\n")
	cat(unique(query_df$Gene_2))
	cat("\n")
}

## Separate the read counts array by sample and use the sampleID as the column name (Sample Metadata)
LecturesArrayToSamples <- function(query_df, CohortMetadata){
	colnames(query_df) <- c("Chrom", "SJ_start", "SJ_end","CohortLectures_array", "Gene_SJ_start", "Gene_SJ_end")
	Samples_index <- CohortMetadata$sample_id

	query_df <- query_df %>% separate(CohortLectures_array, Samples_index, ",") # Annot the samople name by its position in the array
	query_df$SJ_Key <- paste0(query_df$Chrom,"_", query_df$SJ_start, "_", query_df$SJ_end)

	query_df <- query_df %>% select(Chrom, SJ_start, SJ_end, SJ_Key, everything())

	return(query_df)
}

## Format to work with the dataframe
TransposeAnnotator <- function(splitted_lectures_df, CohortMetadata){

	Samples_index <- CohortMetadata$sample_id
	df_t <- as.data.frame(t(splitted_lectures_df[,Samples_index])) # Select columns by sample name (Samples_index) and transpose
	colnames(df_t) <- splitted_lectures_df$SJ_Key # The colums are now the Splice Junctions
	df_t$sample_id <- Samples_index # Select by Sample name

	merged_df <- merge(df_t, CohortMetadata, by.X="sample_id", by.y="sample_id",all.x = TRUE)


	row.names(merged_df) <- merged_df$sample_id #unsorted

	merged_df <- merged_df[order(merged_df$INDEX),] # INDEX number of the Cohort Metatdata csv File (the merge function un-sorts the file)
	merged_df <- merged_df %>% select(INDEX, sample_id, case_id, file_id.BAM, file_name.BAM, "file_id.STAR-SJCounts", "file_name.STAR-SJCounts", everything())

	return(merged_df)
}

## Write query function:
writeQuery <- function(Gene, df, mapping_option, output_directory) {
	if (is.null(mapping_option)) {
		write.table(df, paste0(output_directory, Gene,"_annotSJ.tsv"), sep= "\t", quote=FALSE, row.names=FALSE)
	} else {
		write.table(df, paste0(output_directory, Gene,"_", mapping_option, "_annotSJ.tsv"), sep= "\t", quote=FALSE, row.names=FALSE)
	}
}

## MAIN Wrapping_function
mainQueryProcess <- function(gene, query_df, CohortMetadata, mapping_option, output_directory){
	## Input: splice junction query result
	## Output: sample annotated and transposed df to further work with
	#### 1) Check selected annotations:
	selectedAnnotations(query_df)
	#### 2) Separate the read counts array by sample and use the sampleID as the column name (Sample Metadata)
	splitted_lectures_df <- LecturesArrayToSamples(query_df, CohortMetadata)
	#### 3) Format to work with the dataframe:
	Annot_t_df <- TransposeAnnotator(splitted_lectures_df, CohortMetadata)
	#### 4) Save results
	writeQuery(gene, Annot_t_df, mapping_option, output_directory)
}

if (cohort_project == "BeatAML"){
	## BeatAML
	### Cohort Metadata
	head(CohortMetadata) # The sample's index corresponds to its position in the collected lectures array
	CohortMetadata <- CohortMetadata[order(CohortMetadata$INDEX),] # Ensure the INDEX is sorted
	Samples_index <- CohortMetadata$sample_id # To recover sample information by their index position

	### Splice Junctions files
	cat("Loading splice junction collection....\n")
	AML_SJ <- read.delim(SJ_collection, sep = "\t")

	if (is.null(mapping_option) || mapping_option == "UM"){
		### 1. ASXL1
		Gene="ASXL1"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%ASXL1%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 2. DNMT3A
		Gene="DNMT3A"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%DNMT3A%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 3. EP300
		Gene="EP300"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%EP300%' and Gene_1 not like '%EP300-AS1%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 4. FLT3
		Gene="FLT3"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%FLT3%' and Gene_1 not like '%FLT3LG%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 5. IDH1
		Gene="IDH1"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%IDH1%' and Gene_1 not like '%IDH1-AS1%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 6. KDM6A
		Gene="KDM6A"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%KDM6A%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 7. KMT2D
		Gene="KMT2D"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%KMT2D%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 8. KRAS
		Gene="KRAS"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%KRAS%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 9. NRAS
		Gene="NRAS"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%NRAS%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 10. PTPN11
		Gene="PTPN11"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%PTPN11%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 11.TET2
		Gene="TET2"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%TET2%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 12.TP53
		Gene="TP53"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 = 'TP53(NM_000546);TP53(NM_001126112);TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001126115);TP53(NM_001126116);TP53(NM_001126117);TP53(NM_001126118);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276697);TP53(NM_001276698);TP53(NM_001276699);TP53(NM_001276760);TP53(NM_001276761)' or Gene_1 = 'TP53(NM_000546);TP53(NM_001126112);TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001126118);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276760);TP53(NM_001276761)'or Gene_1 = 'TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001126116);TP53(NM_001126117);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276698);TP53(NM_001276699)'or Gene_1 = 'TP53(NM_000546);TP53(NM_001126112);TP53(NM_001126115);TP53(NM_001126118);TP53(NM_001276697);TP53(NM_001276760);TP53(NM_001276761)'or Gene_1 = 'TP53(NM_001126113);TP53(NM_001126117);TP53(NM_001276695);TP53(NM_001276699)'
		or Gene_1 = 'TP53(NM_001126114);TP53(NM_001126116);TP53(NM_001276696);TP53(NM_001276698)'
		or Gene_1 = 'TP53(NM_000546);TP53(NM_001126112);TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276760);TP53(NM_001276761)' or Gene_1 = 'TP53(NM_001126112);TP53(NM_001276761)' or Gene_1 = 'TP53(NM_000546);TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001126118);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276760)'
		or Gene_1 = 'TP53(NM_000546);TP53(NM_001126112);TP53(NM_001126113);TP53(NM_001126114);TP53(NM_001126118);TP53(NM_001276695);TP53(NM_001276696);TP53(NM_001276760);TP53(NM_001276761);WRAP53(NM_001143990)'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 13. U2AF1
		Gene="U2AF1"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like 'U2AF1(NM_001025203);U2AF1(NM_001025204);U2AF1(NM_006758);U2AF1L5(NM_001320646);U2AF1L5(NM_001320648);U2AF1L5(NM_001320650);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025204);U2AF1(NM_006758);U2AF1L5(NM_001320646);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025203);U2AF1(NM_001025204);U2AF1L5(NM_001320648);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025203);U2AF1L5(NM_001320648)' or Gene_1 like 'U2AF1(NM_001025204);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_006758);U2AF1L5(NM_001320646)'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)

		### 17.WT1
		Gene="WT1"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%WT1%' and Gene_1 not like '%WT1-AS%' and Gene_1 not like '%SWT1%'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)
	}

	if (mapping_option == "MM"){
		Gene="U2AF1"
		cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
		query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like 'U2AF1(NM_001025203);U2AF1(NM_001025204);U2AF1(NM_006758);U2AF1L5(NM_001320646);U2AF1L5(NM_001320648);U2AF1L5(NM_001320650);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025204);U2AF1(NM_006758);U2AF1L5(NM_001320646);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025203);U2AF1(NM_001025204);U2AF1L5(NM_001320648);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_001025203);U2AF1L5(NM_001320648)' or Gene_1 like 'U2AF1(NM_001025204);U2AF1L5(NM_001320651)' or Gene_1 like 'U2AF1(NM_006758);U2AF1L5(NM_001320646)'")
		mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
		rm(Gene,query_df)
	}

}else if(cohort_project =="TCGA"){
	## TCGA
	### Cohort Metadata
	CohortMetadata # The sample's index corresponds to its position in the collected lectures array
	CohortMetadata <- CohortMetadata[order(CohortMetadata$INDEX),] # Ensure the INDEX is sorted
	Samples_index <- CohortMetadata$sample_id # To recover sample information by their index position

	### Splice Junctions files
	cat("Loading splice junction collection....\n")
	AML_SJ <- read.delim(SJ_collection, sep = "\t")

	### 1. NRAS
	Gene="NRAS"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%NRAS%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 2. KRAS
	Gene="KRAS"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%KRAS%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 3. KMT2D
	Gene="KMT2D"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%KMT2D%' and Gene_1 not like '%DDN%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 4. FLT3
	Gene="FLT3"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%FLT3%' and Gene_1 not like '%FLT3LG%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 5. IDH1
	Gene="IDH1"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%IDH1%' and Gene_1 not like '%IDH1-AS1%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 6. WT1
	Gene="WT1"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 like '%WT1%' and Gene_1 not like '%WT1-AS%' and Gene_1 not like '%SWT1%'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)

	### 7. CBL
	Gene="CBL"
	cat(paste0("\n-------------------------------------------------------------------\nExtracting splice junctions of ", Gene, "....\n"))
	query_df <- sqldf("select SJ_1, SJ_2, SJ_3, SJ_4, Gene_1, Gene_2 from AML_SJ where Gene_1 = 'CBL(NM_005188)'")
	mainQueryProcess(Gene, query_df,CohortMetadata, mapping_option, output_directory)
	rm(Gene,query_df)
}else{
	cat("\nProject of interest not defined\n")
}
