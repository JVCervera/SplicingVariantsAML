#!/bin/bash

# Configuration variables

# Cohort variables:

COHORT_PROJECT="TCGA"
INPUT_EXCEL="/home/mguaita/SplicingVariants/Data/TCGA_cohort.xlsx"
VCF_CSV="/home/mguaita/SplicingVariants/Data/TCGA.AMLExomeFiles.csv"
RNASEQ_CSV="/home/mguaita/SplicingVariants/Data/TCGA.AMLRNASeqFiles.csv"
RESULTS_DIR="/home/mguaita/SplicingVariants/Results/TCGA"

# Raw data directories

VCF_DIR="/opt/Hematologia/SplicingVariantsAML/GDC_VCF/TCGA"
BAM_DIR="/path/to/BAM" #Without last "/"
STAR_SJCOUNTS_DIR="/opt/Hematologia/SplicingVariantsAML/GDC_SJ/TCGA"
VARIANTS_TO_SEARCH="/home/mguaita/SplicingVariants/Data/Potential_SAVs.tsv"


# Scripts paths

PYTHON_SCRIPTS_DIR="/home/mguaita/SplicingVariants/Scripts" #Without last "/"
R_SCRIPTS_DIR="/home/mguaita/SplicingVariants/Scripts" #Without last "/"
RNAMUT_EXECUTABLE_DIR="/path/to/rnamut" #Without last "/"
GENE_PANEL="/path/to/rnamut/custom_GenePanel.ind"

# Wxs vcf processing variables

MUTMERGED_OUT="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.mut_merged.txt"
MUTMERGED_ANNOT="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.mut_merged.annot.txt"
VARSEARCH_DNA_RESULTS="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.FoundVariantsDNA.tsv"
VARSEARCH_EXCEL_DNA_VAF_RESULT="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.FoundVariantsDNA.vaf.xlsx"

# Rnaseq variant calling variables

FASTQ_DIR="/path/to/fastq"
RNAMUT_VARCALL_DIR="/home/mguaita/SplicingVariants/Results/TCGA/RNAMut" #Without last "/"
VARSEARCH_EXCEL_RNA_RESULT="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.FoundVariantsRNA.xlsx"

# Rnaseq splice junction processing variables

SJ_COLLECTION_FILE="/home/mguaita/SplicingVariants/Results/TCGA/TCGA.SJcollection"
SJ_OUT_DIR="/home/mguaita/SplicingVariants/Results/TCGA/SpliceJunction"
