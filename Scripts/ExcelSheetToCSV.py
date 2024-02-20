#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to convert a Excel file sheet to a CSV file. Creates the output file in the same folder of the provided excel file.

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

def convert_and_check_columns(input_excel_path, output_csv, sheet_name, header):
	try:
		# Read the Excel file
		df = pd.read_excel(input_excel_path, sheet_name=sheet_name)
		df

		path, filename = os.path.split(input_excel_path)
		csv_output_path = os.path.join(path, output_csv)

		if header == False:
			df.to_csv(csv_output_path, header=False, index=False)
			print("Excel conversion and column check passed successfully.")

		elif header == True:
			df.to_csv(csv_output_path, index=False)
			print("Excel conversion and column check passed successfully.")

	except Exception as e:
		raise ValueError(f"Error converting Excel file: {str(e)}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Script to convert a Excel file sheet to a CSV file. In the same directory of the input excel generates a csv file adapted to force the required column order format: 'sample_id', 'case_id', 'file_id.BAM' ,'file_name.BAM' required for the next processing steps of the pipeline.")
	parser.add_argument('-i', '--inputexcel', required=True, help="Path to the input Excel file.")
	parser.add_argument('-s', '--sheetname', required=True, help="Sheet name of the Excel file to be converted.")
	parser.add_argument('-o', '--outputcsv', required=True, help="Name of the output CSV file.")
	parser.add_argument('-C', '--columnnames', required=False, action='store_true', help="OPTIONAL. Include the header line in the output CSV file.")
	parser.add_argument('-D', '--debug', action='store_true', help=argparse.SUPPRESS)

	args = parser.parse_args()
	# Debug mode:
	if args.debug:
		developer_debug(args)

	try:
		input_excel_path = args.inputexcel
		sheet_name = args.sheetname
		output_csv = args.outputcsv
		if args.columnnames:
			header = args.columnnames

		else:
			header=False

		if args.debug:
			developer_debug(args)

		convert_and_check_columns(input_excel_path, output_csv, sheet_name, header)

	except ValueError as ve:
		print(f"Error: {ve}")
