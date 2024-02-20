#!/bin/bash

# Configuration variables

# Cohort variables:

COHORT_PROJECT="TCGA"
INPUT_EXCEL="/Data/TCGA_cohort.xlsx"
VCF_CSV="/Data/TCGA.AMLExomeFiles.csv"
RNASEQ_CSV="/Data/TCGA.AMLRNASeqFiles.csv"
RESULTS_DIR="/Results/TCGA" #Without last "/"

# Raw data directories

VCF_DIR="/TCGAv36/VCF" #Without last "/"
BAM_DIR="/TCGAv36/STAR2Pass-Genome" #Without last "/"
STAR_SJCOUNTS_DIR="/TCGAv36/STAR-SJCounts" #Without last "/"
VARIANTS_TO_SEARCH="/Data/Potential_SAVs.tsv"


# Scripts paths

PYTHON_SCRIPTS_DIR="/Scripts" #Without last "/"
R_SCRIPTS_DIR="/Scripts" #Without last "/"
RNAMUT_EXECUTABLE_DIR="/path/to/rnamut" #Without last "/"
GENE_PANEL="/path/to/rnamut/custom_GenePanel.ind"

# Wxs vcf processing variables

MUTMERGED_OUT="/Results/TCGA/TCGA.mut_merged.txt"
MUTMERGED_ANNOT="/Results/TCGA/TCGA.mut_merged.annot.txt"
VARSEARCH_DNA_RESULTS="/Results/TCGA/TCGA.FoundVariantsDNA.tsv"
VARSEARCH_EXCEL_DNA_VAF_RESULT="/Results/TCGA/TCGA.FoundVariantsDNA.vaf.xlsx"

# Rnaseq variant calling variables

FASTQ_DIR="/path/to/fastq" #Without last "/"
RNAMUT_VARCALL_DIR="/Results/TCGA/RNAMut" #Without last "/"
VARSEARCH_EXCEL_RNA_RESULT="/Results/TCGA/TCGA.FoundVariantsRNA.xlsx"

# Rnaseq splice junction processing variables

SJ_COLLECTION_FILE="/Results/TCGA/TCGA.SJcollection"
SJ_OUT_DIR="/Results/TCGA/SpliceJunction"
