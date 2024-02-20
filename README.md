# SplicingVariants in AML

## Index

* [Dependencies](#dependencies)
* [Pipeline](#pipeline)
  + [WXS Variant Calling data processing](#WXS-Variant-Calling-data-processing)
  + [RNASeq Variant Calling](#RNASeq-Variant-Calling)
  + [RNASeq Splice Junction Collection](#RNASeq-Splice-Junction-Collection)
  + [Integration](#Integration)
* [Splice altering effect analysis](#Splice-altering-effect-analysis)
* [Contact](#contact)
* [Publication](#publication)


## Dependencies
Required software to run the pipeline:
 + Python 3.0
 + R version 4.1.2
 + bam2fastq version 1.1.0
 + RNAMut version 1.2a
 + junc utils version 0.5.2

### Python packages

+ pandas
+ argparse
+ pysam

### R libraries
+ optparse
+ dplyr
+ tidyr
+ readxl
+ openxlsx
+ stringr
+ sqldf
+ ggplot2
+ DT
+ plotly
+ knitr
+ ggsignif
+ RcmdrMisc

## Bioinformatic Pipeline
The main bash script `Integration_pipeline.sh` incorporates all three workflows to process the different biological data structures.
BeatAML cohort:
  + Configuration file: `Data/BeatAML_config.sh`
  + File listings: `Data/BeatAML_cohort.xlsx` (or `Data/BeatAML.AMLExomeFiles.csv` and `Data/BeatAML.AMLRNASeqFiles.csv`)
```
bash Integration_pipeline.sh BeatAML_config.sh
```

TCGA cohort:
  + Configuration file: `Data/TCGA_config.sh`
  + File listings: `Data/TCGA_cohort.xlsx` (or `Data/TCGA.AMLExomeFiles.csv` and `Data/TCGA.AMLRNASeqFiles.csv`)
```
bash Integration_pipeline.sh TCGA_config.sh 
```
The environment variables and the working directories are sourced from the cohort-specific files, the identifiers of the processed files are listed in the cohort-specific xlsx files (as bash does not support excel, the pipeline starts converting the excel sheets to csv files before running the workflows). The individual scripts can also be run separately. 

###  WXS Variant Calling data Processing 

* **Step 1**. Collection of unique somatic variants

```
usage: python variant_collector_vcf.py [-h]
                                        --inputfile 
                                        --workingdir
                                        --mode 
                                        --outputfile

For additional help: python variant_collector_vcf.py -h    
```

Example of usage with the BeatAML cohort:
```
python variant_collector_vcf.py
-i BeatAML.AMLExomeFiles.csv
-w path/to/BeatAML/VCF
-M4
-o BeatAML.mut_merged.txt
```

* **Step 2**. Sample information recovery
```
usage: Rscript recoverVariant_SampleID.R [-h]
                                          --cohort_project {"BeatAML", "TCGA"}
                                          --cohort_metadata
                                          --variant_collection
                                          --output_file

For additional help: Rscript recoverVariant_SampleID.R -h  
```

Example of usage with the BeatAML cohort:
```
Rscript recoverVariant_SampleID.R 
-p BeatAML
-c BeatAML.AMLExomeFiles.csv
-i BeatAML.mut_merged.txt
-o BeatAML.mut_merged.annot.txt
```

* **Step 3**. Variant searcher in WXS Variant calling data

```
usage: Rscript find_DNAvariants.R [-h]
                                --cohort_project {"BeatAML", "TCGA"}
                                --rna_metadata
                                --search_variant
                                --variant_collection
                                --output_file

For additional help: Rscript find_DNAvariants.R -h  
```

Example of usage with the BeatAML cohort:
```
Rscript find_DNAvariants.R
-p BeatAML
-r BeatAML.AMLRNASeqFiles.csv
-s Potential_SAVs.tsv
-v BeatAML.mut_merged.annot.txt
-o BeatAML.FoundVariantsDNA.tsv
```

* **Step 4**. Get Variant Allele Frequency (VAF)
```
usage: python get_VAF_vcf.py [-h]
                              --inputfile
                              --workingdir 
                              --foundVariants
                              --outputfile
```

Example of usage with the BeatAML cohort:
```
python get_VAF_vcf.py
-i BeatAML
-w path/to/BeatAMLv36/VCF
-f BeatAML.FoundVariantsDNA.tsv
-o BeatAML.FoundVariantsDNA.vaf.xlsx
```

###  RNASeq Variant Calling 

* **Step 5**. BAM to FASTQ conversion

bam2fastq is used for the file conversion.
```
bam2fastq -o  path/to/FASTQ/file_1.fq path/to/FASTQ/file_2.fq path/to/BAM/file.bam
```

* **Step 6**. RNAMut Variant Calling

RNAMut software is used for variant calling.
```
java -jar /path/to/RNAmut.jar
-n sample_id
/path/to/FASTQ/file_1.fq,/path/to/FASTQ/file_2.fq
/path/to/RNAMut_VariantCalling_Results
-i GENE_PANEL
-f oncogenicity_filter.txt 
```

* **Step 7**. Variant searcher in RNASeq Variant calling data
```
usage: python FormatRNAMutOutput.py [-h]
                                     --cohortproject {"BeatAML","TCGA"}
                                     --SAV
                                     --dnametadata
                                     --workingdirectory
                                     --outputfile

For additional help: python FormatRNAMutOutput.py -h
```

Example of usage with the BeatAML cohort:
```
python FormatRNAMutOutput.py 
-c BeatAML
-s Potential_SAVs.tsv
-d BeatAML.AMLExomeFiles.csv
-w path/to/RNAMut_VariantCalling_Results
-o BeatAML.RNAMutResults.xlsx
```

### RNASeq Splice Junction Collection 

* **Step 8**. Splice Junction Collection
```
usage: python SpliceJunctionCollector.py [-h]
                                          --inputfile 
                                          --workingdirectory 
                                          --mappingtype {"UM","MM"} 
                                          --outputfile

For additional help: python SpliceJunctionCollector.py -h
```

Example of usage with the BeatAML cohort:
```
python RNASpliceJunctionCollector.py 
-i BeatAML_cohort.xlsx
-d /STAR_SpliceJunctions
-m "UM" 
-o BeatAML.SJ_UMcollection
```

* **Step 9** Add Gene notation

The gene notation is added to the collected splice junctions using the function annotate of the junc-utils tool.
```
junc_utils annotate --genome_id hg38 BeatAML.SJ_UMcollection.txt BeatAML.SJ_UMcollection_annot.txt
```

* **Step 10** SQL-based queries to extract a particular gene splice junctions

The SQL queries are defined based on the gene notation of the splice junctions
```
usage: Rscript SQL-queries.R [-h]
                              --cohort_project {"BeatAML", "TCGA"}
                              --cohort_metadata
                              --mapping_option
                              --splice_junction
                              --output_directory

For additional help: Rscript SQL-queries.R -h
```

Example of usage with the BeatAML cohort:
```
Rscript SQL-queries.R
-p BeatAML 
-c BeatAML.AMLExomeFiles.csv
-m "UM"
-i BeatAML.SJ_UMcollection
-o /Results/BeatAML/SpliceJunction
```
### Integration
* **Step 11** WXS and RNASeq variant calling results integration
```
usage: Rscript integrate_foundvariants.R [-h]
                                          --dnavariants
                                          --rnavariants
                                          --output

For additional help: Rscript integrate_foundvariants.R -h
```

Example of usage with the BeatAML cohort:
```
Rscript integrate_foundvariants.R
-d BeatAML.FoundVariantsDNA.vaf.xlsx
-r BeatAML.FoundVariantsRNA.xlsx
-o BeatAML.Variant_summary.xlsx
```

## Splice Altering Effect Analysis
The potential splice altering effect analysis is conducted at the `BeatAML_report.Rmd` and  `TCGA_report.Rmd` scripts and both generate an interactive results report.
 + BeatAML results: `Results/BeatAML_report.html`
 + TCGA results: `Results/TCGA_report.html`

## Contact

Jose Cervera 

Hematology Researh Group, Instituto de Investigación Sanitaria La Fe

46026

Valencia, Spain

## Publication (TBD)
