#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to find the variants of interest, in our specific case Papaemmanuil's variants with a potential splice altering effect called (candidate SAV), in the RNAMut variant calling results. Collects all the variants called by RNAMut and if a variant of interest is found is flagged as "CandidateSAV",the search is based on the RNAMut Mutation annotation. Generates two tab-separated output files: 1. All collected variants, 2. Only variants of interest.

"""

# Dependencies and libraries
import os, gzip, subprocess
import pandas as pd
import argparse

def developer_debug(args):
	"""
	Developer debug: print arguments and exits.
	"""
	print("Developer Debug Mode:")
	print("----------------------------------------------------")
	for arg, value in vars(args).items():
		print(f"{arg}: {value}")
	print("Exiting without further execution.")
	print("----------------------------------------------------")
	exit()

# Input cheks:
def check_working_directory(working_directory, substring="mutation-all"):
	try:
		if not os.path.exists(working_directory):
			raise FileNotFoundError(f"Error: Provided RNAMut Directory '{working_directory}' does not exist.")

		if not os.listdir(working_directory):
			raise ValueError(f"Error: Provided RNAMut Directory '{working_directory}' is empty.")

		matching_files = [file for file in os.listdir(working_directory) if substring in file]
		if not matching_files:
			raise ValueError(f"Error: No '{substring}' files found in '{working_directory}'.")

	except Exception as e:
		raise ValueError(str(e))

def create_output_directory(working_directory, dir_to_create):
	# Create somatic folder
	working_directory = working_directory.rstrip('/')
	path = os.path.join(working_directory, dir_to_create)
	os.makedirs(path, exist_ok=True)

def generate_RNAMutKey(input_string):
	# Genomic notation pattern
	pattern = re.compile(r'(\d+):g\.(\d+)([ACGT])>([ACGT])')
	match = re.match(pattern, input_string)

	if match:
		# Matched groups
		chromosome, position, ref_base, alt_base = match.groups()
		# RNAMut_formatHg38
		RNAMutKey = f"chr{chromosome}:{position}-{position}_Sub:{ref_base}>{alt_base}"
		return RNAMutKey
	else:
		# If no match is found, return the original string
		return input_string

def get_RNAMut_variants(df):
	try:
		if 'RNAMut_formatHg38' not in df.columns:
			print("Missing RNAMut notation column\n. Creating it from genomic notation hg38.")
			if 'hg38' not in df.columns:
				raise ValueError(f"The candidate SAV tab-separated file must contain the column {', '.join(required_columns)} with the RNAMut mutation notation or the genomic notation to generate it (column 'hg38')")
			else:
				print("Generating RNAMut notation")
				df['RNAMut_formatHg38'] = df['hg38'].apply(lambda x: generate_RNAMutKey(x))
				variants_list = df["RNAMut_formatHg38"].to_list()
		else:
			variants_list = df["RNAMut_formatHg38"].to_list()

		return variants_list
	except Exception as e:
		raise ValueError(f"Error checking candidate SAV file: {str(e)}")

def load_input_file_data(file_path):
	# Load CSV data and check for RNAMut notation
	SAV = pd.read_csv(file_path, sep='\t')
	SAV_list = get_RNAMut_variants(SAV)
	return SAV_list

##### FUNCTIONS
def get_mutation_files(cohort_project, working_directory):
	list_dir = os.listdir(working_directory)
	mutation_files = [s for s in list_dir if "mutation-all" in s]
	sample_list = []
	pt_list = []
	for s in mutation_files:
		s1 = s.replace("mutation-all_", "")
		s2 = s1.replace(".txt", "")
		sample_list.append(s2)
		if cohort_project == "BeatAML":
			pt_list.append(s2[0:-8]) # remove ".txt" at the end of the filename to get patient ID
		elif cohort_project == "TCGA":
			pt_list.append(s2[0:-4]) # remove "AAA.txt" at the end of the filename to get patient ID

	return mutation_files, sample_list, pt_list

def generate_MutationKey(input_string):
	# Function to convert RNAMut notation to our Mutation_Key
	parts=input_string.split(':')
	chr=parts[0]
	pos=parts[1].split("-")[0]
	ref_alt=parts[2].replace(">",",")

	Mutation_Key=chr+","+pos+","+ref_alt

	return Mutation_Key

def VariantCollector_CandidateSAVearcher(mutation_files, sample_list, pt_list, SAV_list):
	SAMPLE_ID = []
	PATIENT = []
	GENE = [] #F[0]
	MUTATION = [] #F[1]
	PROTMUT = [] #F[2]
	MUTReads = [] #F[3]
	WTReads = [] #F[4]
	VAF = [] #F[5]
	INFO = [] #F[6]
	CandidateSAV = []

	pt_index = 0
	for mt_file in mutation_files:
		with open(mt_file, 'r') as f:
			i=0
			for line in f:
				F = line.rstrip("\n").split("\t")
				if F[0] == "Gene":
					continue
				else:
					SAMPLE_ID.append(sample_list[pt_index])
					PATIENT.append(pt_list[pt_index])
					GENE.append(F[0])
					MUTATION.append(F[1])
					PROTMUT.append(F[2])
					MUTReads.append(F[3])
					WTReads.append(F[4])
					VAF.append(F[5])
					INFO.append(F[6])
					if F[1] in SAV_list:
						CandidateSAV.append("CandidateSAV")
					else:
						CandidateSAV.append("NO")
		pt_index+=1

	df_aux = pd.DataFrame(list(zip(SAMPLE_ID, PATIENT, GENE, MUTATION, PROTMUT, MUTReads, WTReads, VAF, INFO, CandidateSAV)),
					columns = ["SampleID", "Case" ,"Gene", "Mutation", "ProtMut", "MUTReads", "WTReads", "VAF", "Info", "CandidateSAV"])


	return df_aux

def DNA_sample_search(cohort_project,formatted_df, dnaMetadata):
	if cohort_project == "BeatAML":
		dnaMetadata["RNA_Sample"] = dnaMetadata["DNA_Sample"].str.replace("D", "R")

		merged_df = pd.merge(formatted_df, dnaMetadata, left_on='RNA_Sample', right_on='RNA_Sample', how='left')

		merged_df['DNA_Sample'] = merged_df['DNA_Sample'].fillna('No VCF file')

		merged_df.drop(columns=['case_id_y'], inplace=True)
		merged_df = merged_df.rename(columns={'case_id_x': 'case_id'})

	if cohort_project == "TCGA":

		merged_df = pd.merge(formatted_df, dnaMetadata, left_on='case_id', right_on='case_id', how='left')
		merged_df['DNA_Sample'] = merged_df['DNA_Sample'].fillna('No VCF file')

	return merged_df

def save_results(formatted_df, output_file):
	formatted_df.to_excel(output_file, index=False)

def main():

	# ArgumentParser
	parser = argparse.ArgumentParser(description=__doc__)

	parser.add_argument('-c', '--cohortproject', choices=['BeatAML', 'TCGA'], required=True, help='Cohort project defined to set the type of string formating for sample_id: BeatAML (BA0000R) or TCGA (TCGA-AB-0000-00A)')
	parser.add_argument('-s', '--SAV', required=True, help='Path to the candidate SAV tab-separated file')
	parser.add_argument('-d', '--dnametadata', required=True, help='Path to the Cohort Metadata file with the VCF file id and file name identifiers')
	parser.add_argument('-w', '--workingdir', required=True, help='Path to the RNAMut results directory')
	parser.add_argument('-o', '--outputfile', required=True, help='Path and name of the Excel output.')
	parser.add_argument('-D', '--debug', action='store_true', help=argparse.SUPPRESS)

	# Parse the command-line arguments
	args = parser.parse_args()

	# Debug mode:
	if args.debug:
		developer_debug(args)

	cohort_project=args.cohortproject
	input_file = args.SAV
	dna_metadata = args.dnametadata
	working_directory = args.workingdir
	output_file = args.outputfile

	print("----------------------------------------------------")
	print("Cohort project:", cohort_project)
	print("Candidate SAV file:", input_file)
	print("Cohort Metadata CSV file:", dna_metadata)
	print("Path to RNAMut results:", working_directory)
	print("Output Excel file:", output_file)
	print("----------------------------------------------------")

	# RNAMut results files exist:
	working_directory = args.workingdir.rstrip('/') + '/'
	check_working_directory(args.workingdir)
	os.chdir(working_directory)

	# Load SAV data
	try:
		SAV_list = load_input_file_data(input_file)
		print("Loaded candidate SAV data.")
	except ValueError as e:
		print(f"Error: {str(e)}")

	dnaMetadata = pd.read_csv(dna_metadata)
	dnaMetadata = dnaMetadata[['sample_id', 'case_id']]
	dnaMetadata = dnaMetadata.rename(columns={'sample_id': 'DNA_Sample'})

	# RNAMut results format
	mutation_files, sample_list, pt_list = get_mutation_files(cohort_project, working_directory)

	formatted_df = VariantCollector_CandidateSAVearcher(mutation_files, sample_list, pt_list, SAV_list)

	formatted_df['MutationKey_Hg38'] = formatted_df['Mutation'].apply(lambda x: generate_MutationKey(x))
	if cohort_project == "BeatAML":
		formatted_df[['case_id', 'RNA_Sample']] = formatted_df['SampleID'].str.split('_', expand=True)
		formatted_df.drop(columns=["SampleID", "Case"], inplace=True)
		formatted_df['sample_id'] = formatted_df['case_id'] + "_" + formatted_df['RNA_Sample'].str.replace('R', '')

	if cohort_project == "TCGA":
		formatted_df["sample_id"] = formatted_df["Case"]
		formatted_df["case_id"] = formatted_df["Case"]
		formatted_df["RNA_Sample"] = formatted_df["SampleID"]
		formatted_df.drop(columns=["SampleID", "Case"], inplace=True)

	formatted_df["Validable"] = "Validable"
	formatted_df = DNA_sample_search(cohort_project,formatted_df, dnaMetadata)

	formatted_df = formatted_df[["Gene", "MutationKey_Hg38", "sample_id", "case_id", "RNA_Sample", "DNA_Sample", "Validable","Mutation","ProtMut","MUTReads","WTReads","VAF","Info","CandidateSAV"]]

	create_output_directory(working_directory, "Collected_RNAMut")

	out_RNAMut_formattedResults = working_directory + "Collected_RNAMut/RNAMutResults.collection.xlsx"

	save_results(formatted_df, out_RNAMut_formattedResults)

	filtered_df = formatted_df[formatted_df['CandidateSAV'] == 'CandidateSAV']

	save_results(filtered_df, output_file)

	print("Saved results")

if __name__ == "__main__":
	main()
