# find variants of interest
# Script to find a set of variants of interest in a variant collection
# Also takes the GDC cohort manifest to check if the sample harboring the variants presents RNASeq data

required_packages <- c("optparse","dplyr", "readxl", "stringr", "tidyr", "knitr")

# Check if each package is installed and load it if available
#for (package in required_packages) {
#  if (!requireNamespace(package, quietly = TRUE)) {
#    install.packages(package, dependencies = TRUE)  # Install the package if not installed
#  }
#}

  
get_validatable_variants <- function(cohort_project, manifest, foundVariants){
  if (cohort_project == "BeatAML"){
    # Identify the BAcodes
    RNASeq <- c()
    BAM_path <- c()
    Validable <- c()
    N_pt_con_BAM <- c()
    for (i in 1:length(foundVariants$Samples_annot)){
      if (foundVariants$N_samples[i]>1){
        sample_annot <- str_split(foundVariants$Samples_annot[i], ";")[[1]]
        RNASeq_sub <- c()
        BAM_path_sub <- c()
        N_BAM = 0
        for (s in 1:length(sample_annot)){
          sample_string <- str_split(sample_annot[s], "_")[[1]]
          
          caseID <- sample_string[2]
          BAcode <- sample_string[1]
          
          q <- sqldf(paste0("select * from manifest where case_id ='",caseID,"' and sample_id like '%",BAcode,"R%'"))
          q_BAM <- sqldf("select file_id, file_name from q where file_name like '%rna_seq.genomic%'")
          
          if (length(q_BAM$file_id)>=1){
            rna <- "TRUE"
            BAM <- paste0(q_BAM$file_id, "/", q_BAM$file_name)
            N_BAM = N_BAM +1
          }else{
            rna <- "FALSE"
            BAM <- "NO"
          }
          RNASeq_sub <- append(RNASeq_sub, rna)
          BAM_path_sub <- append(BAM_path_sub, BAM)
        }
        RNASeq <- append(RNASeq, paste(RNASeq_sub,collapse=";"))
        BAM_path <- append(BAM_path, paste(BAM_path_sub,collapse=";"))
        N_pt_con_BAM <- append(N_pt_con_BAM, N_BAM)
      }else{
        sample_annot <- foundVariants$Samples_annot[i]
        sample_string <- str_split(sample_annot, "_")[[1]]
        caseID <- sample_string[2]
        BAcode <- sample_string[1]
        N_BAM = 0
        
        q <- sqldf(paste0("select * from manifest where case_id ='",caseID,"' and sample_id like '%",BAcode,"R%'"))
        q_BAM <- sqldf("select file_id, file_name from q where file_name like '%rna_seq.genomic%'")
        
        if (length(q_BAM$file_id)>=1){
          rna <- "TRUE"
          BAM <- paste0(q_BAM$file_id, "/", q_BAM$file_name)
          N_BAM = N_BAM + 1
        }else{
          rna <- "FALSE"
          BAM <- "NO"
        }
        RNASeq <- append(RNASeq,rna)
        BAM_path <- append(BAM_path,BAM)
        N_pt_con_BAM <- append(N_pt_con_BAM, N_BAM)
      }
    } 
    foundVariants$RNASeq <- RNASeq
    foundVariants$BAM <- BAM_path
    foundVariants$N_pt_con_BAM <- N_pt_con_BAM
    
    foundVariants$Validable <- ifelse(grepl("TRUE", foundVariants$RNASeq), "YES", "NO")
    table(foundVariants$Validable == "YES")
  }
  
  if (cohort_project == "TCGA"){
    RNASeq <- c()
    BAM_path <- c()
    Validable <- c()
    N_pt_con_BAM <- c()
    
    for (i in 1:length(foundVariants$Samples_annot)){
      if (foundVariants$N_samples[i]>1){
        sample_annot <- str_split(foundVariants$Samples_annot[i], ";")[[1]]
        RNASeq_sub <- c()
        BAM_path_sub <- c()
        N_BAM = 0
        for (s in 1:length(sample_annot)){
          sample_string <- str_split(sample_annot[s], "_")[[1]]
          
          caseID <- paste((str_split(sample_string[1], "-",4)[[1]][1:3]), collapse="-")
          sample_code <- sample_string[1]
          sample_type <- sample_string[2]
          
          q <- sqldf(paste0("select * from manifest where case_id ='",caseID,"'"))
          q_BAM <- sqldf("select file_id, file_name from q where file_name like '%rna_seq.genomic%'")
          
          if (length(q_BAM$file_id)>=1){
            rna <- "TRUE"
            BAM <- paste0(q_BAM$file_id, "/", q_BAM$file_name)
            N_BAM = N_BAM +1
            
          }else{
            rna <- "FALSE"
            BAM <- "NO"
          }
          RNASeq_sub <- append(RNASeq_sub, rna)
          BAM_path_sub <- append(BAM_path_sub, BAM)
        }
        RNASeq <- append(RNASeq, paste(RNASeq_sub,collapse=";"))
        BAM_path <- append(BAM_path, paste(BAM_path_sub,collapse=";"))
        N_pt_con_BAM <- append(N_pt_con_BAM, N_BAM)
        
      }else{
        sample_annot <- foundVariants$Samples_annot[i]
        sample_string <- str_split(sample_annot, "_")[[1]]
        
        caseID <- paste((str_split(sample_string, "-",4)[[1]][1:3]), collapse="-")
        sample_code <- sample_string[1]
        sample_type <- sample_string[2]
        
        q <- sqldf(paste0("select * from manifest where case_id ='",caseID,"'"))
        q_BAM <- sqldf("select file_id, file_name from q where file_name like '%rna_seq.genomic%'")
        
        N_BAM = 0
        
        if (length(q_BAM$file_id)>=1){
          rna <- "TRUE"
          BAM <- paste0(q_BAM$file_id, "/", q_BAM$file_name)
          N_BAM = N_BAM + 1
          
        }else{
          rna <- "FALSE"
          BAM <- "NO"
        }
        RNASeq <- append(RNASeq,rna)
        BAM_path <- append(BAM_path,BAM)
        N_pt_con_BAM <- append(N_pt_con_BAM, N_BAM)
        
      }
    }  
    foundVariants$RNASeq <- RNASeq
    foundVariants$BAM <- BAM_path
    foundVariants$N_pt_con_BAM <- N_pt_con_BAM
    
    foundVariants$Validable <- ifelse(grepl("TRUE", foundVariants$RNASeq), "YES", "NO")
    table(foundVariants$Validable == "YES")
  }
  return(foundVariants)
}

sample_level_output <- function(cohort_project, foundVariants){
  samples_list <- foundVariants$Samples_annot
  new_rows_list <- list()
  transformed_df <- data.frame(
    Gene = character(),
    Mutation_key = character(),
    sample_id = character(),
    case_id = character(),
    DNA_Sample = character(),
    RNA_Sample = character(),
    Validable = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(samples_list)) {
    for (sample in samples_list[[i]]) {
      # Extract values from Mutation key and Samples
      gene <- foundVariants$SYMBOL[i]
      mutation_key <- foundVariants$Mutation_KeyHg38[i]
      sample_annot <- str_split(foundVariants$Samples_annot[i], ";")[[1]]
      
      for (s in 1:length(sample_annot)){
        if (cohort_project == "BeatAML"){
          sample_string <- str_split(sample_annot[s], "_")[[1]]
          caseID <- sample_string[2]
          BAcode <- gsub("D", "",sample_string[1])
          RNA_sampleid <- manifest$sample_id[manifest$case_id == caseID & manifest$sample_id == paste0(BAcode, "R")]
          if (length(RNA_sampleid) == 0){
            RNA_sampleid = "No RNA BAM"
            Validable_str = "No Validable"
          }else{
            RNA_sampleid = RNA_sampleid
            Validable_str = "Validable"
          }
          
          # Create a new row in the transformed data frame
          new_row <- data.frame(
            Gene = gene,
            Mutation_key = mutation_key,
            sample_id = paste0(caseID, "_", BAcode),
            case_id = caseID,
            DNA_Sample = paste0(BAcode, "D"), # From the mut merged collection , hence DNA sample
            RNA_Sample = RNA_sampleid, # From the rna metadata, hence RNA sample
            Validable = Validable_str,
            stringsAsFactors = FALSE
          )
        }
        
        if (cohort_project == "TCGA"){
          sample_string <- str_split(sample_annot[s], "-")[[1]]
          caseID <- paste(sample_string[1],sample_string[2],sample_string[3],sep="-")
          RNA_sampleid <- manifest$sample_id[manifest$case_id == caseID] 
          if (length(RNA_sampleid) == 0){
            RNA_sampleid = "No RNA BAM"
            Validable_str = "No Validable"
          }else{
            RNA_sampleid = RNA_sampleid
            Validable_str = "Validable"
          }
          # Create a new row in the transformed data frame
          new_row <- data.frame(
            Gene = gene,
            Mutation_key = mutation_key,
            sample_id = caseID,
            case_id = caseID,
            DNA_Sample = sample_annot[s], # From the mut merged collection , hence DNA sample
            RNA_Sample = RNA_sampleid, # From the rna metadata, hence RNA sample
            Validable = Validable_str,
            stringsAsFactors = FALSE
          )
        }
        new_rows_list <- c(new_rows_list, list(new_row))
      }
    }
  }
  
  transformed_df <- do.call(rbind, new_rows_list)
  # Append the new row to the transformed data frame
  return(transformed_df)  
}

save_results <-function(df, output_file){
  write.table(df, file = output_file, sep = "\t", col.names = TRUE, row.names = FALSE)
}

library(optparse)
# Option parsing
option_list <- list(
  make_option(c("-p", "--cohort_project"), type="character", help="Set to BeatAML or TCGA"),
  make_option(c("-r", "--rna_metadata"), type="character", help="Path to cohort metadata csv file with the BAM files identifiers "),
  make_option(c("-s", "--search_variant"), type="character", help="Path to the set of variants to search with hg38 mutation key."),
  make_option(c("-v", "--variant_collection"), type="character", help="Path to the unique somatic variant collection"),
  make_option(c("-o", "--output_file"), type="character", help="Path to the output annotated unique somatic variant collection. Generates a tab-separated file.")
)

opt_parser <- OptionParser(
  description =paste0("R Script to add sample notation to the unique somatic mutation collection"),
  
  usage = "Rscript find_variants.R -p <cohort_project> -r <rna_metadata> -s <search_variant> -v <variant_collection> -o <output_file>",
  
  option_list = option_list)
opt <- parse_args(opt_parser)



# Check if all three arguments are defined
if (is.null(opt$cohort_project) || is.null(opt$rna_metadata) || is.null(opt$search_variant)  || is.null(opt$variant_collection) || is.null(opt$output_file)) {
  cat("Please provide values for all arguments:-p <cohort_project> -r <rna_metadata> -s <search_variant> -v <variant_collection> -o <output_file> \n")
  cat("For help, run: Rscript recoverVariant-sampleID.R -h\n")
  q(save = "no")
}

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

# Get the Cohort Metadata, the variant collection path and the directory path to save the extracted Splice Junctions
cohort_project <- opt$cohort_project
rna_metadata_path <- opt$rna_metadata
search_variant_path <- opt$search_variant
variant_collection_path <- opt$variant_collection
output <- opt$output_file

# Check the Cohort Metadata file
tryCatch({
  rna_metadata <- read_csv(rna_metadata_path)
  manifest <- as.data.frame(rna_metadata) # Not a tibble
  manifest <- manifest[,c("sample_id", "case_id", "file_id.BAM", "file_name.BAM")]
  column_names <- c("sample_id", "case_id", "file_id", "file_name")  
  colnames(manifest) <- column_names 
  
}, error = function(e) {
  stop("Error reading the metadata files.\n")
})


# Load set of variants to search:
tryCatch({
  variants_of_interest <- read.delim(search_variant_path, sep = "\t")
}, error = function(e) {
  stop("Error reading the provided variant set to search file.\n")
})


# Load unique variant collection:
tryCatch({
  mut_merged <- read.delim(variant_collection_path, sep = "\t")
}, error = function(e) {
  stop("Error reading the provided variant collection file.\n")
})

# Find variants:
#table(mut_merged$Mutation_KeyHg38 %in% variants_of_interest$MutationKey_Hg38)[2] # 14 variants in BeatAML / 8 variants TCGA
foundVariants <- mut_merged[mut_merged$Mutation_KeyHg38 %in% variants_of_interest$MutationKey_Hg38,]
foundVariants <- merge(foundVariants, variants_of_interest[, c("MutationKey_Hg38", "SYMBOL")], by.x = "Mutation_KeyHg38", by.y="MutationKey_Hg38")
print(as.data.frame(foundVariants[,c("Mutation_KeyHg38", "SYMBOL")]))

# How many variants can be validated?:
validatable_variants <- get_validatable_variants(cohort_project, manifest, foundVariants) # 9 variants BeatAML / 6 variants TCGA

output_without_ext <- tools::file_path_sans_ext(output)
formatted_output_name <- paste0(output_without_ext, ".variantsummary.")
output_summary <- paste0(formatted_output_name, tools::file_ext(output))

save_results(validatable_variants,output_summary)

foundVariants_formated <- sample_level_output(cohort_project,foundVariants)
save_results(foundVariants_formated,output)
cat(paste0("Found variants saved at: ", output_summary,"\nSample level results saved at: ",output, "\n"))

