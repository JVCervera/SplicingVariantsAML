#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to collect all the unique splice junctions (SJ) from a list of STAR-generated SJ.out.tab files, and extract the uniquely mapping reads or the multi-mapping reads of all the samples.
This script has been adapted from SAVNet software (https://github.com/friend1ws/SAVNet).
The data input is an CSV file, the programme requires the two columns of the GDC Manifest named as "file_id.STAR-SJCounts" and "file_name.STAR-SJCounts" of the selected STAR-SJCounts files.
"""

# Dependencies and libraries
import os, gzip, subprocess, argparse

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

required_libraries = ['pandas', 'argparse']
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

# Input cheks:

def check_csv_columns(file_path):
	try:
		cohortMetadata = pd.read_csv(file_path)
		required_columns = ['file_id.STAR-SJCounts', 'file_name.STAR-SJCounts']

		if not all(column in cohortMetadata.columns for column in required_columns):
			raise ValueError(f"The provided CSV '{file_path}' must contain the columns {', '.join(required_columns)}. Please check column names")

	except Exception as e:
		raise ValueError(f"Error checking CSV file: {str(e)}")

def load_cohortMetadata(file_path):
	try:
	# Load RNASeq file data
		cohortMetadata = pd.read_csv(file_path)
		return cohortMetadata

	except Exception as e:
		raise ValueError(f"Error loading the RNASeq file data from CSV '{file_path}': {str(e)}")


def check_working_directory(cohortMetadata, workingdirectory):
	for index, row in cohortMetadata.iterrows():
		SJ_fileID = row["file_id.STAR-SJCounts"]
		SJ_fileName = row["file_name.STAR-SJCounts"]

		directory_path = os.path.join(workingdirectory, SJ_fileID)
		if not os.path.exists(directory_path):
			raise FileNotFoundError(f"Directory '{directory_path}' does not exist.")

		file_path = os.path.join(directory_path, SJ_fileName)
		if not os.path.exists(file_path):
			raise FileNotFoundError(f"File '{file_path}' does not exist.")

# Splice Junction Functions
#Part of the code adapted from SAVNet:
def extract_unique_SJ(SJ_file_list):
	splice_junctions = {}
	for SJ_file in SJ_file_list:
		hin = gzip.open(SJ_file, "rt")
		for line in hin:
			F = line.rstrip('\n').split('\t')
			if F[0].startswith('#'): continue
			key = F[0] + '\t' + F[1] + '\t' + F[2]
			if key not in splice_junctions: splice_junctions[key] = 1

	return(splice_junctions)

def get_SJ_reads(SJ_file_list,splice_junctions,mappingtype,outputfile):
	temp_id = 0
	hout = open(outputfile + ".tmp.unsorted.txt", 'w')
	for SJ_file in SJ_file_list:
		hin = gzip.open(SJ_file, "rt")
		for line in hin:
			F = line.rstrip('\n').split('\t')
			if F[0].startswith('#'): continue
			if F[0] + '\t' + F[1] + '\t' + F[2] in splice_junctions:
				if mappingtype == "UM":
					print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + str(temp_id) + '\t' + F[6], file = hout) # F[6] - Uniquely mapping; F[7] - Multi-mapping
				elif mappingtype == "MM":
					print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + str(temp_id) + '\t' + F[7], file = hout) # F[6] - Uniquely mapping; F[7] - Multi-mapping

		temp_id = temp_id + 1

	hout.close()

	hout = open(outputfile + '.tmp.sorted.txt', 'w')
	subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", outputfile + ".tmp.unsorted.txt"], stdout = hout)
	hout.close()

	print("Writing Splice Junction Collection....")
	temp_chr = ""
	temp_start = ""
	temp_end = ""
	temp_count = ["0"] * temp_id
	hout = open(outputfile, 'w')
	with open(outputfile + '.tmp.sorted.txt', 'r') as hin:
		for line in hin:
			F = line.rstrip('\n').split('\t')

			if F[1] != temp_start or F[2] != temp_end or F[0] != temp_chr:

				# if not the first line
				if temp_chr != "":
					print(temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count), file = hout)

				temp_chr = F[0]
				temp_start = F[1]
				temp_end = F[2]
				temp_count = ["0"] * temp_id


			temp_count[int(F[3])] = F[4]

	# last check

	if temp_chr != "":
		print(temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count), file = hout)

	hout.close()

	# remove intermediate files
	subprocess.call(["rm", "-rf", outputfile + ".tmp.unsorted.txt"])
	subprocess.call(["rm", "-rf", outputfile + ".tmp.sorted.txt"])


def main():
	# ArgumentParser
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--inputfile', required=True, help='Path to the Cohort Metadata CSV file with the STAR-SJcounts file identifiers')
	parser.add_argument('-d', '--workingdir', required=True, help='Path to STAR-SJCounts directory')
	parser.add_argument('-m', '--mappingtype', choices=['UM', 'MM'], required=True, help='Type of mapping: Uniquely-mapping reads ("UM") or multi-mapping reads ("MM")')
	parser.add_argument('-o', '--outputfile', required=True, help='Path and name of the output splice junction collection: tab-separared txt file')
	parser.add_argument('-D', '--debug', action='store_true', help=argparse.SUPPRESS)

	# Parse the command-line arguments
	args = parser.parse_args()

	if args.debug:
		developer_debug(args)

	# Check for the required sheet and columns
	try:
		check_csv_columns(args.inputfile)
	except ValueError as e:
		print(f"Error: {str(e)}")
		return

	# Set the working directory
	working_directory = args.workingdir.rstrip('/') + '/'
	os.chdir(working_directory)

	input_file = args.inputfile
	mapping_type = args.mappingtype
	output_file = args.outputfile

	print("-----------------------------------------------------")
	print("Cohort Metadata CSV file:", input_file)
	print("Path to STAR-SJCounts directory:", working_directory)
	print("Reads Mapping type:", mapping_type)
	print("Output file:", output_file)
	print("-----------------------------------------------------")


	# Load CSV data
	try:
		cohortMetadata = load_cohortMetadata(input_file)
		print("Loaded CSV data:", cohortMetadata)
	except ValueError as e:
		print(f"Error: {str(e)}")

	# Check for the input files to exist:
	check_working_directory(cohortMetadata, working_directory)

	# SJ Collector
	SJ_fileID = cohortMetadata["file_id.STAR-SJCounts"].tolist()
	SJ_fileName = cohortMetadata["file_name.STAR-SJCounts"].tolist()
	SJ_filelist = [working_directory + SJ_fileID[f] + "/" + SJ_fileName[f] for f in range(0, len(SJ_fileID))]

	print("Collecting Splice Junctions....")
	UniqueSJs = extract_unique_SJ(SJ_filelist)
	print("Collecting supporting reads....")
	get_SJ_reads(SJ_filelist,UniqueSJs,mapping_type,output_file)

if __name__ == "__main__":
	main()
