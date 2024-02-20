#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to extract the somatic variants of a list of VCF files, collect the variants of the four variant callers in a unique file for each sample and
generate a collection of the unique somatic variants of the whole cohort.\n
Three execution modes:\n
	A) VCF file with the somatic variants of a particular variant caller/patient.\n
	B) VCF file with the somatic variants of the four variant callers by patient.\n
	C) VCF file with all the somatic variants of the complete cohort of patients.\n
"""

# Dependencies check
import os, gzip, subprocess

def required_libraries_check(libraries):
	for library in libraries:
		try:
			if library =="pandas":
				import pandas as pd
			else:
				__import__(library)
		except ImportError:
			print(f"The library '{library}' is not installed.")
			try:
				subprocess.call(['pip', 'install', library])
			except subprocess.CalledProcessError:
				raise RuntimeError(f"Failed to install '{library}'. Please install the required libraries manually and run the script again.")

required_libraries = ['pandas', 'argparse', 'pysam']
required_libraries_check(required_libraries)
import pandas as pd
import argparse
import pysam

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
		cohortMetadata = pd.read_csv(file_path)
		return cohortMetadata

	except Exception as e:
		raise ValueError(f"Error loading data from '{file_path}': {str(e)}")

# Functions:
def create_sampleTAG_list(sample_id_list, case_id_list):
	sampleTAGlist = [f"{elem1}_{elem2}" for elem1, elem2 in zip(sample_id_list, case_id_list)]
	return sampleTAGlist

def create_output_directory(working_directory, dir_to_create):
	# Create somatic folder
	working_directory = working_directory.rstrip('/')
	path = os.path.join(working_directory, dir_to_create)
	os.makedirs(path, exist_ok=True)



# 1. Extraction of somatic variants (generate 1 vcf from each patient)
"""Using the relation of VCF files of the AML samples of the cohort, we automate the extraction proccess
of the somatic variants with a specific function for each variant caller (MuTect2, VarScan2, SomaticSniper and MuSE)"""

def get_SOMATIC_PASS_MuTect(sampleID, vcf_file, working_directory):
	# FORMAT OUTPUTS:
	hout_somatic = open(working_directory +"/ExtractedSomatic/" + sampleID + ".MuTect2.somatic.txt", "w")

	## Write headers hout:
	hout_somatic.write("# SOMATIC PASS vcf \n")
	hout_somatic.write("# Variant caller: MuTect2 \n")
	hout_somatic.write("# Inclusion criteria: absence of germline_risk y panel_of_normals flags, FILTER = PASS \n")
	hout_somatic.write("# CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	## VCF file processing:
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "rt")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	for line in hin:
		F = line.rstrip("\n").split("\t")
		if F[0].startswith("#"): continue
		FILTER = F[6]
		# ONLY THE SOMATIC VARIANTS of high-quality: filter PASS
		if (FILTER.find("panel_of_normals") == -1) and (FILTER.find("germline_risk") == -1) and (FILTER.find("PASS") != -1):
			# Si no tienen etiqueta germinal o normal se asume somatica,
			hout_somatic.write(line)

		## The high quality and trust variants are tagged with the "PASS" flag
		##FILTER=<ID=PASS,Description="All filters passed">

		## In this variant caller the somatic variants are not specifically flagged.
		## However, the germinal variants are marked with this flags (we exclude this recorsds):
		##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">

	hin.close()
	hout_somatic.close()


def get_SOMATIC_PASS_SomaticSniper(sampleID, vcf_file, working_directory):

	# FORMAT OUTPUTS:
	hout_somatic = open(working_directory +"/ExtractedSomatic/" + sampleID + ".SomaticSniper.somatic.txt", "w")

	## Write headers hout:
	hout_somatic.write("# SOMATIC PASS vcf \n")
	hout_somatic.write("# Variant caller: SomaticSniper \n")
	hout_somatic.write("# Inclusion criteria: NORMAL SS = 0 (WT), TUMOR SS = 2 (SOMATIC), FILTER = PASS \n")
	hout_somatic.write("# CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


	## VCF Files processing:
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "rt")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")


	for line in hin:
		F = line.rstrip("\n").split("\t")
		if F[0].startswith("#"): continue

		## FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal, 0=wildtype,1=germline,2=somatic,3=LOH,4=unknown">
		## Somatic variants are flagged with:
		## NORMAL SS = 0, TUMOR SS = 2

		## Germinal variants are flagged with:
		## NORMAL SS = 1, TUMOR SS = 2

		NORMAL = F[9]
		NORMAL_l = NORMAL.split(":")
		TUMOR = F[10]
		TUMOR_l = TUMOR.split(":")

		FILTER = F[6]

		if ((NORMAL_l[11] == "0") and (TUMOR_l[11] == "2") and FILTER == "PASS"):
			hout_somatic.write(line)

	hin.close()
	hout_somatic.close()

def get_SOMATIC_PASS_VarScan(sampleID, vcf_file, working_directory):

	# FORMAT OUTPUTS:
	hout_somatic = open(working_directory +"/ExtractedSomatic/" + sampleID + ".VarScan2.somatic.txt", "w")

	## Write headers hout:
	hout_somatic.write("# SOMATIC PASS vcf \n")
	hout_somatic.write("# Variant caller: VarScan2 \n")
	hout_somatic.write("# Inclusion criteria: columna INFO 'SOMATIC' & SS = 2; FILTER = PASS \n")
	hout_somatic.write("# CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


	## VCF Files processing:
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "rt")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	for line in hin:
		F = line.rstrip("\n").split("\t")
		if F[0].startswith("#"): continue

		## Somatic variants are flagged with:
		## INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)">
		## INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">
		INFO = F[7]

		## En la columna de FILTER = PASS
		FILTER = F[6]

		if (INFO.find(";SOMATIC;") != -1) and (INFO.find(";SS=2;") !=-1) and FILTER == "PASS":
			hout_somatic.write(line)

	hin.close()
	hout_somatic.close()

def get_SOMATIC_PASS_MuSE(sampleID, vcf_file, working_directory):
	# FORMAT OUTPUTS:
	hout_somatic = open(working_directory +"/ExtractedSomatic/" + sampleID + ".MuSE.somatic.txt", "w")

	## Write headers hout:
	hout_somatic.write("# SOMATIC PASS vcf \n")
	hout_somatic.write("# Variant caller: MuSE \n")
	hout_somatic.write("# Inclusion criteria: columna INFO 'SOMATIC' & SS = 2; FILTER = PASS \n")
	hout_somatic.write("# CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	## VCF Files processing:
	if vcf_file.endswith(".vcf"):
		hin = open(vcf_file, "rt")
	elif vcf_file.endswith(".vcf.gz"):
		hin = gzip.open(vcf_file, "rt")

	for line in hin:
		F = line.rstrip("\n").split("\t")
		if F[0].startswith("#"): continue

		## Somatic variants are flagged with:
		## FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">
		## INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">
		INFO = F[7]

		TUMOR = F[10]
		TUMOR_l = TUMOR.split(":")

		##FILTER=<ID=PASS,Description="All filters passed">
		##FILTER=<ID=Tier1,Description="Confident level 1">
		##FILTER=<ID=Tier2,Description="Confident level 2">
		##FILTER=<ID=Tier3,Description="Confident level 3">
		##FILTER=<ID=Tier4,Description="Confident level 4">
		##FILTER=<ID=Tier5,Description="Confident level 5">
		FILTER = F[6]

		if (INFO.find("SOMATIC") != -1) and (TUMOR_l[4] == "2") and (FILTER == "PASS"):
			hout_somatic.write(line)

	hin.close()
	hout_somatic.close()

def main_extract_somatic(cohortMetadata, sampleTAGlist, working_directory):
	variant_callers_list = ["MT", "SS","VS", "MuS"]
	for i in range(0, len(sampleTAGlist)):
		sampleTAG = sampleTAGlist[i]
		for vc in variant_callers_list:
			if vc == "MT":
				vcf_file_id = cohortMetadata.loc[i,"file_id.MuTect2"]
				vcf_file_name = cohortMetadata.loc[i,"file_name.MuTect2"]
				vcf_file = working_directory + vcf_file_id + "/" + vcf_file_name
				get_SOMATIC_PASS_MuTect(sampleTAG, vcf_file, working_directory) # The directories have exist already

			if vc == "VS":
				vcf_file_id = cohortMetadata.loc[i,"file_id.VarScan2"]
				vcf_file_name = cohortMetadata.loc[i,"file_name.VarScan2"]
				vcf_file = working_directory + vcf_file_id + "/" + vcf_file_name
				get_SOMATIC_PASS_VarScan(sampleTAG, vcf_file, working_directory)

			if vc == "SS":
				vcf_file_id = cohortMetadata.loc[i,"file_id.SomaticSniper"]
				vcf_file_name = cohortMetadata.loc[i,"file_name.SomaticSniper"]
				vcf_file = working_directory + vcf_file_id + "/" + vcf_file_name
				get_SOMATIC_PASS_SomaticSniper(sampleTAG, vcf_file, working_directory)

			if vc == "MuS":
				vcf_file_id = cohortMetadata.loc[i,"file_id.MuSE"]
				vcf_file_name = cohortMetadata.loc[i,"file_name.MuSE"]
				vcf_file = working_directory + vcf_file_id + "/" + vcf_file_name
				get_SOMATIC_PASS_MuSE(sampleTAG, vcf_file, working_directory)


#2. Collect together all the variants of one patient called by the four variant callers in one file
def mut_annotation(somatic_file, mut_key_list, vc):
	n_rows = 0

	if somatic_file.endswith(".txt"):
		hin = open(somatic_file, "r")
	elif somatic_file.endswith(".txt.gz"):
		hin = gzip.open(somatic_file, "rt")

	included_variants = 0
	for line in hin:
		F = line.rstrip('\n').split('\t')
		if F[0].startswith('#'): continue

		n_rows +=1

		F[2] = "."
		F[5] = "60"
		F[6] = "PASS"
		F[7] = "SOMATIC"

		# Mutation Key constructor:
		mut_key = []
		mut_key.append(F[0]) # CHROM
		mut_key.append(F[1]) # POS
		mut_key.append(F[2]) # ID
		mut_key.append(F[3]) # REF
		mut_key.append(F[4]) # ALT
		mut_key.append(F[5]) # QUAL
		mut_key.append(F[6]) # FILTER
		mut_key.append(F[7]) # INFO

		mut_key_str = "\t".join(mut_key)
		# Include the mutation if it has not been previously found:

		if mut_key_str not in mut_key_list:
			# Start the variant
			included_variants +=1
			mut_key_list.append(mut_key_str)

	hin.close()

	return(mut_key_list)

def write_mergedVCF(out_vcf, mut_key_list, sample_ID):
	## Writer of the merged VCF file of one variant caller
	with open(out_vcf, "w") as hout:

		x = len(mut_key_list)
		y = 0   # Counting to identify the last mutation
		hout.write("# Patient " + sample_ID + " Merged vcf SOMATIC PASS \n")
		hout.write("# Variant callers: MuTect2, VarScan2, SomaticSniper, MuSE \n")
		hout.write("# CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

		for mut in sorted(mut_key_list):
			y += 1
			if y == x:
				hout.write(str(mut))
			else:
				hout.write(str(mut) + "\n")

def main_merge_somatic(sampleTAGlist, somatic_dir, merged_dir):
	variant_callers_list = ["MT", "SS","VS", "MuS"]

	# Loop to unify the VCFs, 1 VCF per sample that collects the unique variants:
	# sample_list - list of samples
	# variant_callers_list - Variant Callers to unify

	## Stats output format: patient | nºVC1 | nºVC2 | nºVC3 | nºVC4 | NºTotal
	# patient --> sample_list

	n=1
	for i in range(0, len(sampleTAGlist)):
		sample = sampleTAGlist[i]
		ms = str(n) + ": " + sample
		n+=1
		mut_key_list = []
		for vc in variant_callers_list:
			if vc  == "MT":
				somatic_file = somatic_dir + "/" + sampleTAGlist[i] + ".MuTect2.somatic.txt"
			if vc  == "VS":
				somatic_file = somatic_dir + "/" + sampleTAGlist[i] + ".VarScan2.somatic.txt"
			if vc == "SS":
				somatic_file = somatic_dir + "/" + sampleTAGlist[i] + ".SomaticSniper.somatic.txt"
			if vc == "MuS":
				somatic_file = somatic_dir + "/" + sampleTAGlist[i] + ".MuSE.somatic.txt"

			mut_key_list = mut_annotation(somatic_file, mut_key_list, vc) # This function overwrites mut_key_list

		merged_somatic_file = merged_dir + "/" + sample + ".merged.txt"
		write_mergedVCF(merged_somatic_file, mut_key_list, sample)

#3. Collection of unique somatic variants of the whole AML samples of the cohort:
def main_unique_variant_collection(mutation_file_list, output_collection):
	## Code extracted and adapted from SAVNet:
	mut2sample = {}
	sample_ind = 0
	for mut_file in mutation_file_list:
		sample_ind = sample_ind + 1
		is_vcf = True if mut_file.endswith(".txt") or mut_file.endswith(".txt.gz") else False
		hin2 = gzip.open(mut_file, 'rt') if mut_file.endswith(".gz") else open(mut_file, 'r')

		for line2 in hin2:
			F2 = line2.rstrip('\n').split('\t')
			if F2[0].startswith('#'): continue
			if F2[0] == "Chr": continue

			if is_vcf == False:
				pos, ref, alt = F2[1], F2[3], F2[4]

				# insertion
				if F2[3] == "-":
					# get the sequence for the reference base
					seq = ""
					for item in pysam.faidx(reference, F2[0] + ":" + str(F2[1]) + "-" + str(F2[1])):
						seq = seq + item.rstrip('\n')
					seq = seq.replace('>', '')
					seq = seq.replace(F2[0] + ":" + str(F2[1]) + "-" + str(F2[1]), '')
					ref, alt = seq, seq + F2[4]

				# deletion
				if F2[4] == "-":
					# get the sequence for the reference base
					seq = ""
					for item in pysam.faidx(reference, F2[0] + ":" + str(int(F2[1]) - 1) + "-" + str(int(F2[1]) - 1)):
						seq = seq + item.rstrip('\n')
					seq = seq.replace('>', '')
					seq = seq.replace(F2[0] + ":" + str(int(F2[1]) - 1) + "-" + str(int(F2[1]) - 1), '')
					pos, ref, alt = str(int(F2[1]) - 1), seq + F2[3], seq

				QUAL = 60
				INFO = "SOMATIC"

				key = '\t'.join([F2[0], pos, '.', ref, alt, str(QUAL), "PASS", INFO])

			else:

				key = '\t'.join(F2[0:8])

			if key not in mut2sample:
				mut2sample[key] = []

			mut2sample[key].append(str(sample_ind))

		hin2.close()

	sample_num = sample_ind

	hout = open(output_collection, 'w')
	for mut in sorted(mut2sample):
		if len(mut2sample[mut]) == sample_num: continue
		print(mut + '\t' + ','.join(mut2sample[mut]), file = hout)

	hout.close()

def main():
	# ArgumentParser
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--inputfile', required=True, help='Path to the Cohort Metadata file with the VCF file id and file name identifiers')
	parser.add_argument('-w', '--workingdir', required=True, help='Path to the VCF files directory')
	parser.add_argument('-o', '--outputfile', help='Path and name of the output variant collection tab-separared txt file')
	parser.add_argument('-D', '--debug', action='store_true', help=argparse.SUPPRESS)
	parser.add_argument('-M', '--mode', type=int, choices=[1, 2, 3, 4], required=True,
					help=('Execution mode. Choose from:\n'
						  '  1. Extract the somatic variants of the individual files\n'
						  '  2. Merge the extracted somatic variants of each sample\n'
						  '  3. Collect the unique somatic variants\n'
						  '  4. All three steps\n'
						  'If execution mode 3 or 4 please specify the output file path and name (-o).'
						  ))


	# Parse the command-line arguments
	args = parser.parse_args()

	if args.debug:
		developer_debug(args)

	if args.mode in [3, 4]:
		if not args.outputfile:
			parser.error("Please specify output file.")


	# Set the working directory
	working_directory = args.workingdir.rstrip('/') + '/'
	check_working_directory(working_directory)
	os.chdir(working_directory)

	# Value access by attribute names
	input_file = args.inputfile
	mode = args.mode
	if args.outputfile:
		output_collection = args.outputfile

	print("-----------------------------------------------------")
	print("Cohort Metadata CSV file: ", input_file)
	print("Path to VCF directory: ", working_directory)
	print("Execution mode: ", mode)
	if args.outputfile:
		print("Output file: ", output_collection)
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
	sample_id_list = cohortMetadata["sample_id"].to_list()
	case_id_list = cohortMetadata["case_id"].to_list()
	sampleTAGlist = create_sampleTAG_list(sample_id_list, case_id_list)

	# 1. Extract somatic variantas of all the files:
	if mode == "1" or "4":
		print("Extracting the somatic variants of individual files ...")
		create_output_directory(working_directory, "ExtractedSomatic")
		main_extract_somatic(cohortMetadata, sampleTAGlist, working_directory)

	# 2. Merge somatic variants of each individual sample:
	if mode == "2" or "4":
		print("Merging the somatic variants of individual samples ...")
		create_output_directory(working_directory, "MergedSomatic")
		somatic_dir = working_directory + "ExtractedSomatic" # without last "/"
		merged_dir = working_directory + "MergedSomatic" # without last "/"
		main_merge_somatic(sampleTAGlist, somatic_dir, merged_dir)

	# 3. Generate somatic variant collection
	if mode == "3" or "4":
		print("Generating unique somatic variant collection ...")
		create_output_directory(working_directory, "VariantCollection")
		merged_dir = working_directory + "MergedSomatic" # without last "/"
		mutation_file_list = [merged_dir + "/" + sample + ".merged.txt" for sample in sampleTAGlist]
		main_unique_variant_collection(mutation_file_list, output_collection)


if __name__ == "__main__":
	main()
