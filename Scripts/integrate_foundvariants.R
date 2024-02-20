# Script Description
# Script to merge the results of the found variants of interest in DNA and RNA variant calling data.

library(optparse)
# Option parsing
option_list <- list(
  make_option(c("-d", "--dnavariants"), type="character", default=NULL, help="Path to the excel file with the DNA results."),
  make_option(c("-r", "--rnavariants"), type="character", default=NULL, help="Path to the excel file with the RNA results."),
  make_option(c("-o", "--output"), type="character", help="Path to the excel output file.")
)

opt_parser <- OptionParser(
  description =paste0("R Script to integrate the found variants in DNA data with the found variants in RNA data"),
  
  usage = "Rscript Integrate_foundvariants.R -d <dna_variants> -r <rna_variants> -o <output>",
  
  option_list = option_list)
opt <- parse_args(opt_parser)



# Check if all three arguments are defined
if (is.null(opt$cohort_project) || is.null(opt$cohort_metadata) || is.null(opt$splice_junction) || is.null(opt$output_directory)) {
	cat("Please provide values for all four arguments: -d <dna_variants> -r <rna_variants> -o <output> .\n")
	cat("For help, run: Rscript Integrate_foundvariants.R -h\n")
	q(save = "no")
}

# Set working environment:
## Load libraries
library(dplyr)
library(readr)
library(stringr)
library(readxl)
library(openxlsx)
library(tidyr)
library(knitr)

#sessionInfo()

# Get the Cohort Metadata, the SJ collection path and the directory path to save the extracted Splice Junctions
dna_variants_path  <- opt$dnavariants
rna_variants_path  <- opt$rnavariants
output_file   <- opt$output

# Read inputs
dna_variants <- read_excel(dna_variants_path, sheet = "Sheet1")
rna_variants <- read_excel(rna_variants_path, sheet = "Sheet1")

# Same format
dna_variants$CalledBy <- "DNA"
rna_variants$CalledBy <- "RNA"

# Integrate
m <- merge(dna_variants, rna_variants, by=c("MutationKey_Hg38", "Gene", "sample_id", "case_id", "RNA_Sample", "DNA_Sample"), suffixes=c(".DNA", ".RNA"),all= TRUE)
m$CalledBy <- ifelse(m$CalledBy.DNA =="DNA" & is.na(m$CalledBy.RNA), "DNA",
                     ifelse(m$CalledBy.RNA =="RNA" & is.na(m$CalledBy.DNA), "RNA", 
                     ifelse(m$CalledBy.DNA =="DNA" &m$CalledBy.RNA =="RNA", "DNA/RNA",
                     "None")))
m$CalledBy.DNA <- NULL
m$CalledBy.RNA <- NULL
m$Info <- NULL
m$CandidateSAV <- NULL

# Annotate
for(n in 1:(nrow(m))){

  print(n)
  if (m$DNA_Sample[n] == "No VCF file"){
    print(n)
    m$MuTect2_vaf[n] <- m$DNA_Sample[n] #"No VCF File"
    m$VarScan2_vaf[n] <- m$DNA_Sample[n]
    m$SomaticSniper_vaf[n] <- m$DNA_Sample[n]
    m$MuSE_vaf[n] <- m$DNA_Sample[n]
    
    m$MuTect2_filter[n] <- m$DNA_Sample[n] #"No VCF File"
    m$VarScan2_filter[n] <- m$DNA_Sample[n]
    m$SomaticSniper_filter[n] <- m$DNA_Sample[n]
    m$MuSE_filter[n] <- m$DNA_Sample[n]
  }else{# There is VCF DNA Sample
    if (is.na(m$MuTect2_vaf[n])){
      m$MuTect2_vaf[n] <- "Not DNA Called" #"No VCF File"
      m$MuTect2_filter[n] <- "Not DNA Called" #"No VCF File"
    }
    if (is.na(m$VarScan2_vaf[n])){
      m$VarScan2_vaf[n] <- "Not DNA Called" #"No VCF File"
      m$VarScan2_filter[n] <- "Not DNA Called" #"No VCF File"
    }
    if (is.na(m$SomaticSniper_vaf[n])){
      m$SomaticSniper_vaf[n] <- "Not DNA Called" #"No VCF File"
      m$SomaticSniper_filter[n] <- "Not DNA Called" #"No VCF File"
    }
    if (is.na(m$MuSE_vaf[n])){
      m$MuSE_vaf[n] <- "Not DNA Called" #"No VCF File"
      m$MuSE_filter[n] <- "Not DNA Called" #"No VCF File"
    }
  }
}

for(n in 1:(nrow(m))){
  m$RNA_Sample[n]
  if (m$RNA_Sample[n] == "No RNA BAM"){
    m$Mutation[n] <- m$RNA_Sample[n]
    m$ProtMut[n] <- m$RNA_Sample[n]
    m$MUTReads[n] <- m$RNA_Sample[n]
    m$WTReads[n] <- m$RNA_Sample[n]
    m$VAF[n] <- m$RNA_Sample[n]
  }else{
    if(is.na(m$Mutation[n])){
      m$Mutation[n] <- "Not RNA Called"
      m$ProtMut[n] <- "Not RNA Called"
      m$MUTReads[n] <- "Not RNA Called"
      m$WTReads[n] <- "Not RNA Called"
      m$VAF[n] <- "Not RNA Called"    
    }
  }
  
}
m$Validable <- ifelse(m$RNA_Sample == "No RNA BAM", "No Validable", "Validable")
  
m$Validable.DNA <- NULL
m$Validable.RNA <- NULL

write.xlsx(m, output, rowNames=FALSE)

