#!/bin/bash

# Check configuration file is provided
if [ $# -eq 0 ]; then
	echo "Usage: $0 <config_file>"
	echo "Please provide config.sh"
	exit 1
fi

# Get the configuration file path
config_file="$1"

# Check if the config file exists
if [ ! -f "$config_file" ]; then
	echo "Error: Configuration file '$config_file' not found."
	exit 1
fi


echo "$(date) $(hostname)" >&2

## Load pipeline variables
source "$config_file"
python3 --version 
R --version

echo "Formating VCF Files metadata from Excel.... "
python3 $PYTHON_SCRIPTS_DIR/ExcelSheetToCSV.py \
-i $INPUT_EXCEL \
-s "AML Exome Files" \
-o $VCF_CSV \
-C

echo "Formating RNASeq Files Metadata from Excel...."
python3 $PYTHON_SCRIPTS_DIR/ExcelSheetToCSV.py \
-i $INPUT_EXCEL \
-s "AML-RNASeq-Files" \
-o $RNASEQ_CSV \
-C

### WXS VCF Processing
##############################################################################################################################################################################
### Bash script
### Important: Modify paths and input configuration files as necessary
### Set the cohort project
### All input vcf files must be in the same directory
### 1. Extract all the unique somatic variants of the cohort
### 2. Recover sample information and annotate samples harboring the collected variants, format to one allele change per row
### 3. Search for the variants of interest
##############################################################################################################################################################################
echo "Running Variant VCF Collector...."

python3 $PYTHON_SCRIPTS_DIR/variant_collector_vcf.py -i $VCF_CSV -w $VCF_DIR -M4 -o $MUTMERGED_OUT

echo "Annotating samples and formatting DNA variant collection...."
Rscript $R_SCRIPTS_DIR/recoverVariant_sampleID.R \
-p $COHORT_PROJECT \
-c $VCF_CSV \
-i $MUTMERGED_OUT \
-o $MUTMERGED_ANNOT

echo "Potential variants search...."
Rscript find_DNAvariants.R \
-p $COHORT_PROJECT \
-r $RNASEQ_CSV \
-s $VARIANTS_TO_SEARCH \
-v $MUTMERGED_ANNOT \
-o $VARSEARCH_DNA_RESULTS

echo "Get DNA vaf"

python3 $PYTHON_SCRIPTS_DIR/get_DNA_VAF.py -i $VCF_CSV -w $VCF_DIR -f $VARSEARCH_DNA_RESULTS -o $VARSEARCH_EXCEL_DNA_VAF_RESULT

#### RNA Variant Calling
##############################################################################################################################################################################
### Bash script
### Important: Modify paths and input configuration files as necessary
### Set the cohort project
### As bash does not support Excel, the first step is to convert the input cohort configuration excel to a comma-separated csv file
### Required headerles csv file with  column order for the script to work: sample_id, case_id, file_id.BAM, file_name.BAM
### Set the variants to be searched
###
### 1. Convert AML-RNASeq-Files sheet of the Cohort excel file to csv
### 2. BAM to FASTQ conversion
### 3. Run RNAMut software
### 3. RNAMut results processing and find candidate SAV
##############################################################################################################################################################################

echo "BAM to FASTQ conversion and RNAMut Variant Calling...."

sed '1d' $RNASEQ_CSV > $RNASEQ_CSV.noheader.temp

readarray -t sampleID < <(cut -d',' -f1 $RNASEQ_CSV.noheader.temp)
readarray -t caseID < <(cut -d',' -f2 $RNASEQ_CSV.noheader.temp)
readarray -t sampleBAMSid < <(cut -d',' -f3 $RNASEQ_CSV.noheader.temp)
readarray -t sampleBAMSname < <(cut -d',' -f4 $RNASEQ_CSV.noheader.temp)

rm $RNASEQ_CSV.noheader.temp

for i in $(seq 0 $((${#sampleID[@]} - 1)))
do
	sample_id=${sampleID[$i]}
	case_id=${caseID[$i]}
	sampleTAG=${case_id}_${sample_id}

	BAM_RNA_id=${sampleBAMSid[$i]}
	BAM_RNA_name=${sampleBAMSname[$i]}
	BAM_RNA=${BAM_RNA_id}/${BAM_RNA_name}
	fileName=$(basename "$BAM_RNA_name" ".bam")

	echo $i
	echo ${sampleTAG}
	echo ${BAM_RNA}
	echo ${fileName}
		
	# 2. BAM 2 FASTQ

	bam2fastq -o $FASTQ_DIR/${fileName}"#.fq" $BAM_DIR/${BAM_RNA}
	
	# 3. RUN RNAMut on Custom Panel

	java -jar $RNAMUT_EXECUTABLE_DIR/RNAmut.jar -n ${sampleTAG} $FASTQ_DIR/${fileName}_1.fq,$FASTQ_DIR/${fileName}_2.fq $RNAMUT_VARCALL_DIR -i $GENE_PANEL -f $RNAMUT_EXECUTABLE_DIR/oncogenicity_filter.txt
	rm ${fileName}_1.fq ${fileName}_2.fq ${fileName}_M.fq
done

# 4. RNAMut results processing and find candidate SAV
echo "Potential variants search...."
python $PYTHON_SCRIPTS_DIR/FormatRNAMutOutput.py -c $COHORT_PROJECT -s $VARIANTS_TO_SEARCH -d $VCF_CSV -w $RNAMUT_VARCALL_DIR -o $VARSEARCH_EXCEL_RNA_RESULT

#### RNA Splice Junction Processing
##############################################################################################################################################################################
## 1. Collect the unique splice junctions of the cohort and extract their supporting reads of each sample
## 2. Junction annotation by Junc Utils (reference Hg38 genome)
## 3. SQL-based queries to extract splice junctions of the genes harboring the variants of interest
##############################################################################################################################################################################
# Run the RNASpliceJunctionCollector.py script
echo "Extracting uniquely mapping reads"

python3 $PYTHON_SCRIPTS_DIR/SpliceJunctionCollector.py -i $RNASEQ_CSV -d $STAR_SJCOUNTS_DIR -m UM -o ${SJ_COLLECTION_FILE}_UM.txt

# Run Junc utils
echo "Annotating SJ with Junc utils"

junc_utils annotate --genome_id hg38 ${SJ_COLLECTION_FILE}_UM.txt ${SJ_COLLECTION_FILE}_UM.annot.txt

# Splice Junction selection
echo "Extracting Splice Junctions of interest" # Each query has been tailored for each gene's annotation

Rscript $R_SCRIPTS_DIR/SQL-Queries.R -p $COHORT_PROJECT -c $RNASEQ_CSV -m "UM" -i ${SJ_COLLECTION_FILE}_UM.annot.txt -o $SJ_OUT_DIR

if [ "$COHORT_PROJECT" == "BeatAML" ]; then
	echo "Extracting Multi-mapping reads..."

	python $PYTHON_SCRIPTS_DIR/SpliceJunctionCollector.py -i $RNASEQ_CSV -d $STAR_SJCOUNTS_DIR -m MM -o ${SJ_COLLECTION_FILE}_MM.txt

	# Run Junc utils
	echo "Annotating SJ with Junc utils"

	junc_utils annotate --genome_id hg38 ${SJ_COLLECTION_FILE}_MM.txt ${SJ_COLLECTION_FILE}_MM.annot.txt

	# Splice Junction selection
	echo "Extracting Splice Junctions of interest"

	Rscript $R_SCRIPTS_DIR/SQL-Queries.R -p $COHORT_PROJECT -c $RNASEQ_CSV -m "MM" -i ${SJ_COLLECTION_FILE}_MM.annot.txt -o $SJ_OUT_DIR

fi

#### Found Variants Integration
Rscript $R_SCRIPTS_DIR/integrate_foundvariants.R -d $VARSEARCH_EXCEL_DNA_VAF_RESULT  -r $VARSEARCH_EXCEL_RNA_RESULT  -o $RESULTS_DIR/$COHORT_PROJECT.Variant_summary.xlsx

echo "Pipeline Completed"
echo "$(date) $(hostname)" >&2

