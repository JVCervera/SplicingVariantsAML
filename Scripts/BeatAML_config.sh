#!/bin/bash

# Configuration variables

# Cohort variables:

COHORT_PROJECT="BeatAML"
INPUT_EXCEL="/Data/BeatAML_cohort.xlsx"
VCF_CSV="/Data/BeatAML.AMLExomeFiles.csv"
RNASEQ_CSV="/Data/BeatAML.AMLRNASeqFiles.csv"
RESULTS_DIR="/Results/BeatAML" #Without last "/"

# Raw data directories

VCF_DIR="/BeatAMLv36/VCF" #Without last "/"
BAM_DIR="/BeatAMLv36/STAR2pass_Genome" #Without last "/"
STAR_SJCOUNTS_DIR="/BeatAMLv36/STAR-SJCounts" #Without last "/"
VARIANTS_TO_SEARCH="/Data/Potential_SAVs.tsv"


# Scripts paths

PYTHON_SCRIPTS_DIR="/Scripts" #Without last "/"
R_SCRIPTS_DIR="/Scripts" #Without last "/"
RNAMUT_EXECUTABLE_DIR="/path/to/rnamut" #Without last "/"
GENE_PANEL="/path/to/rnamut/custom_GenePanel.ind"


# Wxs vcf processing variables

MUTMERGED_OUT="/Results/BeatAML/BeatAML.mut_merged.txt"
MUTMERGED_ANNOT="/Results/BeatAML/BeatAML.mut_merged.annot.txt"
VARSEARCH_DNA_RESULTS="/Results/BeatAML/BeatAML.FoundVariantsDNA.tsv"
VARSEARCH_EXCEL_DNA_VAF_RESULT="/Results/BeatAML/BeatAML.FoundVariantsDNA.vaf.xlsx"

# Rnaseq variant calling variables

FASTQ_DIR="/path/to/fastq" #Without last "/"
RNAMUT_VARCALL_DIR="/Results/BeatAML/RNAMut" #Without last "/"
VARSEARCH_EXCEL_RNA_RESULT="/Results/BeatAML/BeatAML.FoundVariantsRNA.xlsx"

# Rnaseq splice junction processing variables

SJ_COLLECTION_FILE="/Results/BeatAML/BeatAML.SJcollection"
SJ_OUT_DIR="/Results/BeatAML/SpliceJunction"
