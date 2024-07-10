#### Librerias
library(readxl)
library(stringr)

#### Funciones
getPSI <- function(events.df, SJCounts_dir, outPSI_dir){
SplicingKey = events.df$SplicingKey
SplicingClass = events.df$SplicingClass
GeneSJfile = events.df$GeneSJfile
ExclusionSJ  = events.df$ExclusionSJ
MutationKeyHg38 <- events.df$MutationKeyHg38

samples_v <- c()
PSIvalues_v <- c()
for (n in 1:length(GeneSJfile)){
  print(n)
  GeneSJ = read_excel(paste0(SJCounts_dir, "/", GeneSJfile[n]), sheet= "Sheet 1")
  samples = paste0(GeneSJ$sample_id,collapse = ";")
  samples_v <- append(samples_v, samples)
  
  ExclusionSJ = str_replace(str_replace(events.df$ExclusionSJ[n], ":","_"),"-", "_")
  ExclusionSJ_reads <- GeneSJ[,grepl(ExclusionSJ, colnames(GeneSJ))][[1]]
  
  if (SplicingClass[n] == "Exon Skipping"){
    InclusionSJ_1 <- str_replace(str_replace(str_split_fixed(string = events.df$InclusionSJ[n], pattern = "\\|", 2)[1],":", "_"),"-", "_")
    InclusionSJ_2 <- str_replace(str_replace(str_split_fixed(string = events.df$InclusionSJ[n], pattern = "\\|", 2)[2],":", "_"),"-", "_")
    
    Inclusion_reads <- GeneSJ[,grepl(InclusionSJ_1, colnames(GeneSJ))][[1]] + GeneSJ[,grepl(InclusionSJ_2, colnames(GeneSJ))][[1]]
    
    # The inlucison is the reference, thus: 
    GeneSJ$PSI = round((Inclusion_reads / (ExclusionSJ_reads + Inclusion_reads)),4)
    
  }else if(SplicingClass[n] == "Double Exon Skipping"){
    InclusionSJ_1 <- str_replace(str_replace(str_split_fixed(string = events.df$InclusionSJ[n], pattern = "\\|", 3)[1],":", "_"),"-", "_")
    InclusionSJ_2 <- str_replace(str_replace(str_split_fixed(string = events.df$InclusionSJ[n], pattern = "\\|", 3)[2],":", "_"),"-", "_")
    InclusionSJ_3 <- str_replace(str_replace(str_split_fixed(string = events.df$InclusionSJ[n], pattern = "\\|", 3)[3],":", "_"),"-", "_")
    
    Inclusion_reads <- GeneSJ[,grepl(InclusionSJ_1, colnames(GeneSJ))][[1]] + GeneSJ[,grepl(InclusionSJ_2, colnames(GeneSJ))][[1]] +  GeneSJ[,grepl(InclusionSJ_3, colnames(GeneSJ))][[1]]
    
    # The inlucison is the reference, thus: 
    GeneSJ$PSI = round((Inclusion_reads / (ExclusionSJ_reads + Inclusion_reads)),4)
  
  }else{
    InclusionSJ_1 <- str_replace(str_replace(events.df$InclusionSJ[n], ":","_"),"-", "_")
    Inclusion_reads <- GeneSJ[,grepl(InclusionSJ_1, colnames(GeneSJ))][[1]] 
    
    GeneSJ$PSI = round(Inclusion_reads / (ExclusionSJ_reads + Inclusion_reads),4)
  }
  
  PSIvalues_v <- append(PSIvalues_v, paste(GeneSJ$PSI, collapse = ";"))
  writexl::write_xlsx(GeneSJ, paste0(outPSI_dir, "/",  str_replace_all(MutationKeyHg38[n], ",", "-"), "_",str_replace(SplicingKey[n],":","-") , ".xlsx"))
}  

events.df$samples <- samples_v
events.df$PSIvalues <- PSIvalues_v

return(events.df)
}

meanPSI_WT <- function(events.cohort, outPSI_dir){
  meanPSI_WT_v <- c()
  for (n in 1:nrow(events.cohort)){
    print(n)
    GeneSJ = read_excel(paste0(outPSI_dir, "/", str_replace_all(events.cohort$MutationKeyHg38[n], ",","-"),"_",str_replace(events.cohort$SplicingKey[n],":","-"), ".xlsx"))
    meanPSI_WT <- round(mean(GeneSJ$PSI[GeneSJ$GROUP == "WT"], na.rm=TRUE),4)
    meanPSI_WT_v <- append(meanPSI_WT_v, meanPSI_WT)
  }
  
  return(meanPSI_WT_v)
}


getPSI_MUT <- function(events.cohort, outPSI_dir,cohort){
  MUT_samples_v  <- c()
  PSIvalues_MUT_v <- c()
  for (n in 1:nrow(events.cohort)){
    print(n)
    GeneSJ = read_excel(paste0(outPSI_dir, "/", str_replace_all(events.cohort$MutationKeyHg38[n], ",","-"),"_",str_replace(events.cohort$SplicingKey[n],":","-"), ".xlsx"))
    if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,C"){
      mut_status = "R132G"
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
      
    }else if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,A"){
      mut_status = "R132C"
      MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
    }else if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,T"){
      mut_status = "R132S"
      MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
    }else if (events.df$MutationKeyHg38[n] =="chr11,119278165,G,C"){
      MUTsamples = "TCGA-AB-2914"
    }else if (events.df$MutationKeyHg38[n] =="chr11,119278164,A,T"){
      MUTsamples = "TCGA-AB-2956"
    }else{
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$GROUP == "MUT"] 
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$GROUP == "MUT"] 
      }
      
    }
    
    MUTsamples_s = paste0(MUTsamples,collapse = ";")
    MUT_samples_v <- append(MUT_samples_v, MUTsamples_s)
    if (cohort == "TCGA"){
      PSIvalues_MUT <- paste0(GeneSJ$PSI[GeneSJ$case_id %in% MUTsamples], collapse = ";")  
    }else{
      PSIvalues_MUT <- paste0(GeneSJ$PSI[GeneSJ$sample_id %in% MUTsamples], collapse = ";")  
    }
    
    PSIvalues_MUT_v <- append(PSIvalues_MUT_v, PSIvalues_MUT)
  }
  return(list(PSIvalues_MUT_v, MUT_samples_v))
}

getdeltaPSI_MUT <- function(events.cohort, outPSI_dir, cohort){
  MUT_samples_v  <- c()
  deltaPSIvalues_MUT_v <- c()
  for (n in 1:nrow(events.cohort)){
    print(n)
    GeneSJ = read_excel(paste0(outPSI_dir, "/", str_replace_all(events.cohort$MutationKeyHg38[n], ",","-"),"_",str_replace(events.cohort$SplicingKey[n],":","-"), ".xlsx"))
    meanPSI_WT <- events.cohort$meanPSI_WT[n]
    if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,C"){
      mut_status = "R132G"
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
    }else if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,A"){
      mut_status = "R132C"
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
    }else if (events.df$MutationKeyHg38[n] =="chr2,208248389,G,T"){
      mut_status = "R132S"
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$Variant_status == mut_status]  
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$Variant_status == mut_status]  
      }
    }else if (events.df$MutationKeyHg38[n] =="chr11,119278165,G,C"){
      MUTsamples = "TCGA-AB-2914"
    }else if (events.df$MutationKeyHg38[n] =="chr11,119278164,A,T"){
      MUTsamples = "TCGA-AB-2956"
    }else{
      if (cohort == "TCGA"){
        MUTsamples = GeneSJ$case_id[GeneSJ$GROUP == "MUT"] 
      }else{
        MUTsamples = GeneSJ$sample_id[GeneSJ$GROUP == "MUT"] 
      }
    }
    
    MUTsamples_s = paste0(MUTsamples,collapse = ";")
    MUT_samples_v <- append(MUT_samples_v, MUTsamples_s)
    
    if (cohort == "TCGA"){
      deltaPSIvalues_MUT <- paste0(round((GeneSJ$PSI[GeneSJ$case_id %in% MUTsamples] - meanPSI_WT),4), collapse = ";")
    }else{
      deltaPSIvalues_MUT <- paste0(round((GeneSJ$PSI[GeneSJ$sample_id %in% MUTsamples] - meanPSI_WT),4), collapse = ";")  
    }
    deltaPSIvalues_MUT_v <- append(deltaPSIvalues_MUT_v, deltaPSIvalues_MUT)
  }
  
  return(deltaPSIvalues_MUT_v)
}


# TCGA manual PSI
events.df <- read_excel("/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/PSIconfig.xlsx")
events.df <- events.df[events.df$Cohort =="TCGA",c("Gene","MutationKeyHg38","SplicingKey","SplicingClass","InclusionSJ","ExclusionSJ","GeneSJfile")]
SJCounts_dir = "/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/TCGA/SpliceJunction/NormalizedExpression"
outPSI_dir = "/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/TCGA/SpliceJunction/PSI"
events.df <- getPSI(events.df, SJCounts_dir,outPSI_dir)
#writexl::write_xlsx(events.df, paste0(outPSI_dir,"/TCGA.psi2.xlsx"))

events.tcga <- events.df

events.tcga$meanPSI_WT <- meanPSI_WT(events.tcga,outPSI_dir)
events.tcga[,c("PSI_MUT", "MUTSamples")] <- getPSI_MUT(events.tcga,outPSI_dir, "TCGA")
events.tcga$deltaPSI_MUT <- getdeltaPSI_MUT(events.tcga,outPSI_dir, "TCGA")
col_order <- c("Gene", "MutationKeyHg38","SplicingKey","SplicingClass","InclusionSJ","ExclusionSJ","GeneSJfile","samples","PSIvalues","meanPSI_WT","PSI_MUT","deltaPSI_MUT","MUTSamples")
events.tcga <- events.tcga[,col_order]
#writexl::write_xlsx(events.tcga, paste0(outPSI_dir,"/TCGA.psi2.xlsx"))

# Beat manual PSI
events.df <- read_excel("/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/PSIconfig.xlsx")
events.df <- events.df[events.df$Cohort =="BeatAML",c("Gene","MutationKeyHg38","SplicingKey","SplicingClass","InclusionSJ","ExclusionSJ","GeneSJfile")]
SJCounts_dir = "/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/BeatAML/SpliceJunction/NormalizedExpression"
outPSI_dir = "/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/BeatAML/SpliceJunction/PSI"
events.df <- getPSI(events.df, SJCounts_dir,outPSI_dir)
#writexl::write_xlsx(events.df, paste0(outPSI_dir,"/BeatAML.psi2.xlsx"))

events.beat <- events.df

events.beat$meanPSI_WT <- meanPSI_WT(events.beat,outPSI_dir)
events.beat[,c("PSI_MUT", "MUTSamples")] <- getPSI_MUT(events.beat,outPSI_dir, "BeatAML")
events.beat$deltaPSI_MUT <- getdeltaPSI_MUT(events.beat,outPSI_dir, "BeatAML")
col_order <- c("Gene", "MutationKeyHg38","SplicingKey","SplicingClass","InclusionSJ","ExclusionSJ","GeneSJfile","samples","PSIvalues","meanPSI_WT","PSI_MUT","deltaPSI_MUT","MUTSamples")
events.beat <- events.beat[,col_order]
#writexl::write_xlsx(events.beat, paste0(outPSI_dir,"/BeatAML.psi2.xlsx"))



########## get deltaPSI
library(dplyr)
library(readxl)
library(writexl)

dirpath <- "/media/adminiis/HematoLaFe/SplicingVariantsAML-main/Results/BeatAML/SpliceJunction/PSI"

files <- list.files(dirpath)
chr_files <- files[grepl("^chr", files)]

for (fi in chr_files) {
  
  filepath <- file.path(dirpath, fi)
  df <- read_excel(filepath)
  
  mean_WT_PSI <- mean(df$PSI[df$GROUP == "WT"], na.rm=TRUE)
  
  df$deltaPSI <- df$PSI - mean_WT_PSI
  write_xlsx(df, path = filepath)
}
