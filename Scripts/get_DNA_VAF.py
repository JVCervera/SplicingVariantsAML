#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to get the VAF values for the found variants of interest
"""

# Dependencies and libraries
import os, gzip, subprocess
import pandas as pd
import argparse
from pandasql import sqldf
pysqldf = lambda q: sqldf(q, globals())
from more_itertools import locate

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


# Inputs check
def input_format_check(input_path):
	required_columns = ['sample_id', 'case_id', "file_id.MuTect2", "file_name.MuTect2","file_id.VarScan2","file_name.VarScan2","file_id.SomaticSniper","file_name.SomaticSniper","file_id.MuSE","file_name.MuSE"]

	try:
		cohortMetadata = pd.read_csv(input_path)

		if not all(column in cohortMetadata.columns for column in required_columns):
			missing_columns = [column for column in required_columns if column not in cohortMetadata.columns]

			if missing_columns:
				print(f"Error: The following columns are missing: {', '.join(missing_columns)}")
			else:
				for index, row in cohortMetadata.iterrows():
					print(row)

	except FileNotFoundError:
		print("Error: The specified input file was not found.")

	except pd.errors.EmptyDataError:
		print("Error: The provided file is empty.")

	except Exception as e:
		print(f"An error occurred: {e}")

def check_working_directory(working_dir):
	try:
		if os.path.exists(working_dir) and os.path.isdir(working_dir):
			if not os.listdir(working_dir):
				print(f"The provided working directory at '{working_dir}' is empty.")
		else:
			print(f"The provided working directory at '{working_dir}' does not exist.")
	except Exception as e:
		print(f"An error occurred: {e}")

def load_cohortMetadata(file_path):
	try:
		Metadata = pd.read_csv(file_path)
		return Metadata

	except Exception as e:
		raise ValueError(f"Error loading data from '{file_path}': {str(e)}")

def load_foundVariants(found_file):
	try:
		fvariants = pd.read_csv(found_file, sep='\t')
		return fvariants

	except Exception as e:
		raise ValueError(f"Error loading data from '{file_path}': {str(e)}")

# Functions
def create_sampleTAG_list(sample_id_list, case_id_list):
	sampleTAGlist = [f"{elem1}_{elem2}" for elem1, elem2 in zip(sample_id_list, case_id_list)]
	return sampleTAGlist

def find_indices(list_to_check, item_to_find):
	indices = locate(list_to_check, lambda x: x == item_to_find)
	return list(indices)

def get_VAF_MuTect2(variant, vcf_file):
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "r")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	foundvariant_bool = False
	for line in hin:
		F = line.rstrip('\n').split('\t')
		if F[0].startswith('#'): continue
		# Mutation Key constructor:
		mut_key = []
		mut_key.append(F[0]) # CHROM
		mut_key.append(F[1]) # POS
		mut_key.append(F[3]) # REF
		mut_key.append(F[4]) # ALT
		mut_key_str = ",".join(mut_key)

		if mut_key_str == variant:
			vcf_filter = F[6]
			vaf = F[10].split(":")[2] # Column tumor: AF value
			foundvariant_bool = True

	if foundvariant_bool == False:
		vcf_filter="Not DNA called"
		vaf="Not DNA called"

	return vaf, vcf_filter

def get_VAF_VarScan2(variant, vcf_file):
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "r")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	foundvariant_bool = False
	for line in hin:
		F = line.rstrip('\n').split('\t')
		if F[0].startswith('#'): continue
		# Mutation Key constructor:
		mut_key = []
		mut_key.append(F[0]) # CHROM
		mut_key.append(F[1]) # POS
		mut_key.append(F[3]) # REF
		mut_key.append(F[4]) # ALT
		mut_key_str = ",".join(mut_key)

		if mut_key_str == variant:
			vcf_filter = F[6]
			vaf_percent = F[10].split(":")[-2] #Column tumor: FREQ value percentage <"Variant allele frequency">
			vaf =  float(vaf_percent.rstrip('%')) / 100.0 # return decimal value
			foundvariant_bool = True

	if foundvariant_bool == False:
		vcf_filter="Not DNA called"
		vaf="Not DNA called"

	return vaf, vcf_filter

def get_VAF_SomaticSniper(variant, vcf_file):
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "r")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	foundvariant_bool = False
	for line in hin:
		F = line.rstrip('\n').split('\t')
		if F[0].startswith('#'): continue
		# Mutation Key constructor:
		mut_key = []
		mut_key.append(F[0]) # CHROM
		mut_key.append(F[1]) # POS
		mut_key.append(F[3]) # REF
		mut_key.append(F[4]) # ALT
		mut_key_str = ",".join(mut_key)

		if mut_key_str == variant:
			vcf_filter = F[6]
			vaf_str = F[10].split(":")[3].split(",") # Column tumor: DP4 values <"high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
			vaf_int = [int(n) for n in vaf_str]
			vaf = round((vaf_int[2] + vaf_int[3])/sum(vaf_int), 4) #(alt-forward + alt-reverse bases)/ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
			foundvariant_bool = True

	if foundvariant_bool == False:
		vcf_filter="Not DNA called"
		vaf="Not DNA called"

	return vaf, vcf_filter

def get_VAF_MuSE(variant, vcf_file):
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "r")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	foundvariant_bool = False
	for line in hin:
		F = line.rstrip('\n').split('\t')
		if F[0].startswith('#'): continue
		# Mutation Key constructor:
		mut_key = []
		mut_key.append(F[0]) # CHROM
		mut_key.append(F[1]) # POS
		mut_key.append(F[3]) # REF
		mut_key.append(F[4]) # ALT
		mut_key_str = ",".join(mut_key)

		if mut_key_str == variant:
			vcf_filter = F[6]
			vaf_str = F[10].split(":")[1].split(",") # Column tumor: AD <"Allelic depths for the ref and alt alleles in the order listed">
			vaf_int = [int(n) for n in vaf_str]
			vaf = round((vaf_int[1]/sum(vaf_int)), 4)
			foundvariant_bool = True

	if foundvariant_bool == False:
		vcf_filter="Not DNA called"
		vaf="Not DNA called"

	return vaf, vcf_filter

def main_get_VAF(variant, vcf_file, vc):
	# vaf can be a number or "Not 'DNA' called"
	if vc == "MuTect2":
		vaf, vcf_filter = get_VAF_MuTect2(variant, vcf_file)
	if vc == "VarScan2":
		vaf, vcf_filter = get_VAF_VarScan2(variant, vcf_file)
	if vc == "SomaticSniper":
		vaf, vcf_filter = get_VAF_SomaticSniper(variant, vcf_file)
	if vc == "MuSE":
		vaf, vcf_filter = get_VAF_MuSE(variant, vcf_file)

	return vaf, vcf_filter

def main():
	# ArgumentParser
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--inputfile', required=True, help='Path to the Cohort Metadata file with the VCF file id and file name identifiers')
	parser.add_argument('-w', '--workingdir', required=True, help='Path to the VCF files directory')
	parser.add_argument('-f', '--foundVariants', required=True, help='Path to the found Variants tab separated files')
	parser.add_argument('-o', '--outputfile', help='Path and name of the output vaf Excel file')
	parser.add_argument('-D', '--debug', action='store_true', help=argparse.SUPPRESS)

	# Parse the command-line arguments
	args = parser.parse_args()

	if args.debug:
		developer_debug(args)

	# Set the working directory
	working_directory = args.workingdir.rstrip('/') + '/'
	check_working_directory(working_directory)
	os.chdir(working_directory)

	# Access by attribute names
	input_file = args.inputfile
	found_file = args.foundVariants
	output_file = args.outputfile

	print("-----------------------------------------------------")
	print("Cohort Metadata CSV file: ", input_file)
	print("Path to VCF directory: ", working_directory)
	print("Found variants TSV file: ", found_file)
	print("Output Excel file: ", output_file)
	print("-----------------------------------------------------")

	# Check for the required sheet and columns
	try:
		input_format_check(args.inputfile)
	except ValueError as e:
		print(f"Error: {str(e)}")
		return

	# Load Cohort data
	try:
		cohortMetadata = load_cohortMetadata(input_file)
		print("Loaded Cohort data:", cohortMetadata)
	except ValueError as e:
		print(f"Error: {str(e)}")

	cohortMetadata = load_cohortMetadata(input_file)

	foundVariants = load_foundVariants(found_file)
	print(foundVariants.head())
	gene_list = foundVariants["Gene"].to_list()
	DNA_sample_list = foundVariants["DNA_Sample"].to_list()
	case_id_list = foundVariants["case_id"].to_list()
	sampleTAGlist = create_sampleTAG_list(DNA_sample_list, case_id_list)
	variants_list = foundVariants["Mutation_key"].to_list()

	varcallers = ["MuTect2", "VarScan2", "SomaticSniper", "MuSE"]

	MT_vaf = []
	MT_filter = []
	VS_vaf = []
	VS_filter = []
	SS_vaf = []
	SS_filter = []
	MS_vaf = []
	MS_filter = []

	for v in range(0, len(variants_list)):
		variant = variants_list[v]
		gene = gene_list[v]
		sample_id = DNA_sample_list[v]
		sampleMetadata = cohortMetadata[cohortMetadata['sample_id'] == sample_id]
		for vc in varcallers:
			column_id = "file_id."+ vc
			column_name = "file_name." + vc
			vcf_file_id  = cohortMetadata.loc[cohortMetadata["sample_id"] == sample_id, column_id].values[0]
			vcf_file_name = cohortMetadata.loc[cohortMetadata["sample_id"] == sample_id, column_name].values[0]
			vcf_file = working_directory + vcf_file_id + "/" + vcf_file_name
			vaf, vcf_filter = main_get_VAF(variant, vcf_file, vc)

			if vc == "MuTect2":
				MT_vaf.append(vaf)
				MT_filter.append(vcf_filter)
			if vc == "VarScan2":
				VS_vaf.append(vaf)
				VS_filter.append(vcf_filter)
			if vc == "SomaticSniper":
				SS_vaf.append(vaf)
				SS_filter.append(vcf_filter)
			if vc == "MuSE":
				MS_vaf.append(vaf)
				MS_filter.append(vcf_filter)

	# Store results 
	foundVariants['MuTect2_vaf'] = MT_vaf
	foundVariants['VarScan2_vaf'] = VS_vaf
	foundVariants['SomaticSniper_vaf'] = SS_vaf
	foundVariants['MuSE_vaf'] = MS_vaf

	foundVariants['MuTect2_filter'] = MT_filter
	foundVariants['VarScan2_filter'] = VS_filter
	foundVariants['SomaticSniper_filter'] = SS_filter
	foundVariants['MuSE_filter'] = MS_filter

	foundVariants = foundVariants.rename(columns={'Mutation_key': 'MutationKey_Hg38'})

	# Write excel results
	foundVariants.to_excel(output_file, index=False)

if __name__ == "__main__":
	main()
