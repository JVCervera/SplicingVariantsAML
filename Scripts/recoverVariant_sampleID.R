# recoverVariant-SampleID.Rmd
# Script to add sample notation to the unique somatic mutation collection.
# The collection VCF the samples are identified by an Index number, with this script we change it to the cohort's ID
# List of required packages

required_packages <- c("optparse","dplyr", "sqldf", "readxl", "stringr", "tidyr", "knitr")

# Check if each package is installed and load it if available
#for (package in required_packages) {
#  if (!requireNamespace(package, quietly = TRUE)) {
#    install.packages(package)  # Install the package if not installed
#  }
#}


library(optparse)
# Option parsing
option_list <- list(
  make_option(c("-p", "--cohort_project"), type="character", help="Set to BeatAML or TCGA"),
  make_option(c("-c", "--cohort_metadata"), type="character", help="Path to cohort metadata csv file with the VCF files identifiers "),
  make_option(c("-i", "--variant_collection"), type="character", help="Path to the unique somatic variant collection"),
  make_option(c("-o", "--output_file"), type="character", help="Path to the output annotated unique somatic variant collection. Generates a tab-separated file. Needed for comma-separated mutation Key")
)

opt_parser <- OptionParser(
  description =paste0("R Script to add sample notation to the unique somatic mutation collection"),

  usage = "Rscript recoverVariant-sampleID.R -p <cohort_project>  -c <cohort_metadata> -i <splice_junction> -o <output_file>",

  option_list = option_list)
opt <- parse_args(opt_parser)




# Set working environment:
## Load libraries
library(dplyr)
library(readr)
library(stringr)
library(sqldf)
library(tidyr)
library(knitr)
library(data.table)


#sessionInfo()

# Check if all three arguments are defined
if (is.null(opt$cohort_project) || is.null(opt$cohort_metadata) || is.null(opt$variant_collection) || is.null(opt$output_file)) {
  cat("Please provide values for all four arguments: -p <cohort_project>  -c <cohort_metadata> -i <variant_collection> -o <output_file>.\n")
  if (is.null(opt$cohort_project)){
    cat("Provide cohort project.\n")
  }
  if (is.null(opt$cohort_metadata)){
    cat("Provide cohort_metadata. \n")
  }
  if (is.null(opt$variant_collection)){
    cat("Provide variant_collection. n")
  }
  if (is.null(opt$output_file)){
    cat("Provide output_file\n")
  }
  cat("For help, run: Rscript recoverVariant-sampleID.R -h\n")
  q(save = "no")
}

# Get the arguments
cohort_project <- opt$cohort_project
cohort_metadata <- opt$cohort_metadata
variant_collection_path <- opt$variant_collection
output_file <- opt$output_file

# Check the Cohort Metadata file
tryCatch({
  CohortMetadata <- read_csv(cohort_metadata)
  CohortMetadata <- as.data.frame(CohortMetadata) # Not a tibble
  # Check if the required columns exist
  if ("sample_id" %in% colnames(CohortMetadata) && "INDEX" %in% colnames(CohortMetadata) && "case_id" %in% colnames(CohortMetadata)) {
      cat("The metadata file has been successfully read, and it contains the required sheet and columns.\n")
  } else {
      stop("The metadata file does not contain the required columns (sample_id and/or INDEX).\n")
    }
  }, error = function(e) {
    stop("Error reading the metadata file. Make sure it is a valid metadata file with the specified sample_id, case_id and INDEX columns.\n")
  })

# Load unique variant collection:
tryCatch({
  mut_merged <- read.delim(variant_collection_path, comment.char = "#", sep = "\t", header=FALSE, col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "SAMPLES"))
  }, error = function(e) {
    stop("Error reading the provided variant collection file.\n")
  })

# Processing Step:
# 1. Sample annotation
BAcodes <- c()
Patients <- c()
N_samples <- c()
N_patients <- c()

MM_samples <- mut_merged$SAMPLES

for (c in 1:length(MM_samples)){
  samples_string <- MM_samples[c]
  samples_list <- str_split(samples_string, ",")[[1]]
  N_samples <- append(N_samples, length(samples_list))

  BA_string <- c()
  pt_string <- c()
  N_pt = 0

  if (cohort_project =="BeatAML"){
    for (s in 1:length(samples_list)){
      index = as.integer(samples_list[s])
      sample_id = CohortMetadata$sample_id[CohortMetadata$INDEX == index]
      pt <-  CohortMetadata$case_id[CohortMetadata$INDEX == index]
      BA <- paste0(sample_id, "_", pt)
      Tissue <- CohortMetadata$tissue[CohortMetadata$INDEX == index]
      Stage <- CohortMetadata$stage[CohortMetadata$INDEX == index]
      
      BA_stringtag <- paste0(BA,"_",Stage,"_",Tissue)
      BA_string <- append(BA_string, BA_stringtag)

      if(!pt %in% pt_string){
        pt_string <- c(pt_string, pt)
        N_pt = N_pt +1
      }else{
          pt_string<- pt_string}
    }
  }

  if (cohort_project == "TCGA"){
    for (s in 1:length(samples_list)){
      index = samples_list[s]
      BA <- CohortMetadata$sample_id[CohortMetadata$INDEX == index]
      pt <-  CohortMetadata$case_id[CohortMetadata$INDEX == index]

      BA_string <- append(BA_string, BA)

      if(!pt %in% pt_string){
        pt_string <- c(pt_string, pt)
        N_pt = N_pt +1
      }else{
        pt_string<- pt_string}
    }
  }
    BA_string <- paste(BA_string, collapse=";")
  pt_string <- paste(pt_string, collapse = ";")


  BAcodes <- append(BAcodes, BA_string)
  Patients <- append(Patients, pt_string)
  N_patients <- append(N_patients, N_pt)
}

mut_merged$Samples_annot <- BAcodes
mut_merged$N_samples <- N_samples
mut_merged$Patients <- Patients
mut_merged$N_patients <- N_patients

mut_merged$Mut_KeyHg38 <- paste0(mut_merged$CHROM,",", mut_merged$POS, ",", mut_merged$REF, ",", mut_merged$ALT)

output_without_ext <- tools::file_path_sans_ext(output_file)
formatted_output_name <- paste0(output_without_ext, ".temp.")
output_temp<- paste0(formatted_output_name, tools::file_ext(output_file))

write.table(mut_merged, output_temp, sep= "\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

# 2. Variant annotation formatting for entries with more than one alternative allele (1nt change per row):
## One registry per nucleotide change:
q <- sqldf("select * from mut_merged where ALT like '%,%'")

q1 <- q
q2 <- q

q1$ALT <- gsub(",.*", "", q1$ALT) # select ALT allele 1
q2$ALT <- gsub(".*,","", q2$ALT) # select ALT allele 2

q_grep <- mut_merged[-grep(",",mut_merged$ALT), ] # remove the entries with double ALT allele

mut_merged_1nt <- rbind(q_grep, setnames(q1, names(q_grep)), setnames(q2, names(q_grep)))
mut_merged_1nt$Mutation_KeyHg38 <- paste0(mut_merged_1nt$CHROM,",", mut_merged_1nt$POS,",", mut_merged_1nt$REF, ",", mut_merged_1nt$ALT)

#length(unique(mut_merged_1nt$Mut_KeyHg38)) # Duplicated Mutation Key
#length(unique(mut_merged_1nt$Mutation_KeyHg38)) # New Mutation Key

dupmut_merged <- mut_merged_1nt$Mutation_KeyHg38[duplicated(mut_merged_1nt$Mutation_KeyHg38)]

mut_merged_1nt$DUPLICATED_VARIANT <- ifelse(mut_merged_1nt$Mutation_KeyHg38 %in% dupmut_merged, "TRUE", "FALSE")

dups <- mut_merged_1nt[mut_merged_1nt$DUPLICATED_VARIANT == "TRUE",]
counts <- sqldf("select Mutation_KeyHg38, count(*) as cnt from dups group by Mutation_KeyHg38")
#length(unique(counts$Mutation_KeyHg38))
#table(counts$cnt == "1")
#table(counts$cnt == "2")
#table(counts$cnt == "3")
#table(counts$cnt == "4")

### BeatAML 173492 variants  - (794*1) - (47*2) - (1*3)


#### Data cleaning step: Unify the duplicated information of the registries
df_groupby <- sqldf("select CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, group_concat(SAMPLES) SAMPLES, group_concat(Samples_annot) Samples_annot, sum(N_samples) N_samples, group_concat(Patients) Patients, sum(N_patients) N_patients, Mutation_KeyHg38, DUPLICATED_VARIANT from dups group by Mutation_KeyHg38")

remove_duplicates <- function(row){
  row <- str_replace(row, ";", ",")
  row <- str_split(row, ",")[[1]]
  row <- unique(row)
  row <- paste(row, collapse=",")
  return (row)
}

len <- function(cell){
  cell <- str_split(cell, ",")[[1]]
  n <- length(cell)
  return(n)
}

df_groupby$SAMPLES  <- lapply(df_groupby$SAMPLES, remove_duplicates)
df_groupby$Samples_annot <- lapply(df_groupby$Samples_annot, remove_duplicates)
df_groupby$Patients <- lapply(df_groupby$Patients, remove_duplicates)
df_groupby$N_samples <-  lapply(df_groupby$SAMPLES, len)
df_groupby$N_patients <- lapply(df_groupby$Patients, len)
df_groupby <- apply(df_groupby,2, as.character) 
head(df_groupby)

#### Remove the duplicated entries from the original dataframe
mut_merged_1nt <- mut_merged_1nt[mut_merged_1nt$DUPLICATED_VARIANT == "FALSE",]
mut_merged_1nt$Mut_KeyHg38 <- NULL

#### Incorporate the new registries:
mut_merged_1nt_corrected <- rbind(mut_merged_1nt, df_groupby)
mut_merged_1nt_corrected$DUPLICATED_VARIANT <- NULL  # Remove duplication flag

write.table(mut_merged_1nt_corrected, output_file,
            sep = "\t", quote = FALSE, row.names =  FALSE)


