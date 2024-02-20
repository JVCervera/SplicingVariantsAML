#!/bin/bash

# Configuration variables

# Cohort variables:

COHORT_PROJECT="BeatAML"
INPUT_EXCEL="/home/mguaita/SplicingVariants/Data/BeatAML_cohort.xlsx"
VCF_CSV="/home/mguaita/SplicingVariants/Data/BeatAML.AMLExomeFiles.csv"
RNASEQ_CSV="/home/mguaita/SplicingVariants/Data/BeatAML.AMLRNASeqFiles.csv"
RESULTS_DIR="/home/mguaita/SplicingVariants/Results/BeatAML"

# Raw data directories

VCF_DIR="/opt/Hematologia/SplicingVariantsAML/GDC_VCF/BeatAML"
BAM_DIR="/opt/Hematologia/BeatAMLv36/STAR2pass_Genome" #Without last "/"
STAR_SJCOUNTS_DIR="/opt/Hematologia/SplicingVariantsAML/GDC_SJ/BeatAML"
VARIANTS_TO_SEARCH="/home/mguaita/SplicingVariants/Data/Potential_SAVs.tsv"


# Scripts paths

PYTHON_SCRIPTS_DIR="/home/mguaita/SplicingVariants/Scripts" #Without last "/"
R_SCRIPTS_DIR="/home/mguaita/SplicingVariants/Scripts" #Without last "/"
RNAMUT_EXECUTABLE_DIR="/path/to/rnamut" #Without last "/"
GENE_PANEL="/path/to/rnamut/custom_GenePanel.ind"


# Wxs vcf processing variables

MUTMERGED_OUT="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.mut_merged.txt"
MUTMERGED_ANNOT="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.mut_merged.annot.txt"
VARSEARCH_DNA_RESULTS="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.FoundVariantsDNA.tsv"
VARSEARCH_EXCEL_DNA_VAF_RESULT="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.FoundVariantsDNA.vaf.xlsx"

# Rnaseq variant calling variables

FASTQ_DIR="/path/to/fastq"
RNAMUT_VARCALL_DIR="/home/mguaita/SplicingVariants/Results/BeatAML/RNAMut" #Without last "/"
VARSEARCH_EXCEL_RNA_RESULT="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.FoundVariantsRNA.xlsx"

# Rnaseq splice junction processing variables

SJ_COLLECTION_FILE="/home/mguaita/SplicingVariants/Results/BeatAML/BeatAML.SJcollection"
SJ_OUT_DIR="/home/mguaita/SplicingVariants/Results/BeatAML/SpliceJunction"
