#!/usr/bin/python

# Script to process hgWiggle output (for one region) and calculate PhastCons values (working with Python 2.7.5).
# Note that this script could be iterated together with hgWiggle using a bash script to process full BED files.
# Required modules = re and sys.
# Usage : python Supplemental_File_S4.py -dataset name -mode 1 -wiggle hgWiggle_output
#	dataset = prefix for the output files.
#	mode = select the option to calculate conservation values, the value could be 1 (full region and binned region) or 2 (full region only).
#	wiggle = hgWiggle output file.

# ------------------------------------
# python modules
# ------------------------------------

import re
import sys

# ------------------------------------
# Misc functions
# ------------------------------------

# Function to separate data from a string.
def process_line(f_line):
	for i in range(len(f_line)):
		if f_line[i].isdigit():
			continue
		else:
			val1 = f_line[0:i]
			val2 = f_line[i+1:len(f_line)]
			break
	return val1, val2	

# Function to read the input file, extract values and sort them accordingly to the size of the region.
def process_file(file_name):
	wig_file = open(file_name)
	for line in wig_file:
		# Get the chromosome and positions for the region of interest.
		if line.startswith("#\tchr"):
			chrom = line[19:len(line) - 1]
		elif line.startswith("#\tposition"): 
			pos = line[22:len(line) - 1]
			positions = process_line(pos)
			# Iniatialize matrix for the positions and values.
			size = int(positions[1]) - int(positions[0]) + 1
			phastcons = [0 for x in range(size)]
		elif line[0].isdigit():
			values = process_line(line)
			# Fill the matrix.
			key = int(values[0])
			phastcons[key - int(positions[0]) + 1] = float(values[1])
	wig_file.close()
	data = {'chrom': chrom, 'start': positions[0], 'end': positions[1], 'size': size, 'phastcons': phastcons}
	return data

# ------------------------------------
# Main function
# ------------------------------------
# Read parameters
argv = sys.argv[1:]
while argv:
	argument = argv.pop(0)
	if(re.match('^-dataset', argument)):
		file_name = argv.pop(0)
	elif(re.match('^-mode', argument)):
		mode = argv.pop(0)
	elif(re.match('^-wiggle', argument)):
		wigfile = argv.pop(0)

# Check if the input file contains data and process it or write zero values to represent no conservation.
with open(wigfile, "r") as wig_file:
	if wig_file.read():
		region_data = process_file(wigfile)
	else:
		region_data = {'chrom': 'chrUn', 'start': '1000', 'end': '9999', 'size': 0, 'phastcons': [0]}

wig_file.close()

name = region_data['chrom'] + ":" + region_data['start'] + "-" + region_data['end'] + "\t"

# Select between the two different modes for the output.
# 1) Average conservation of the region and conservation calculated in 50 bins.
# 2) Only the average conservation of the whole region.
if int(mode) == 1:
	# Calculate the bin size.
	bin_size = region_data['size'] / 50
	# Calculate the average phastcons for each bin and the average for all the region.
	avg_phastcons = [0 for x in range(50)]
	n = 0
	i = 0
	if len(region_data['phastcons']) > 1:
		while n < 50:
			bin_values = region_data['phastcons'][i:i + bin_size]
			avg_phastcons[n] = sum(bin_values) / len(bin_values)
			n += 1
			i = i + bin_size
	avg_region = sum(region_data['phastcons']) / len(region_data['phastcons'])
	# Write the results in two different files.
	# File with average conservation for the whole region.
	with open(file_name + "_summary_region.txt", "a+") as sum_file:
		if sum_file.read():
			sum_file.write("\n" + name + str(avg_region))
		else:
			sum_file.write(name + str(avg_region))
	sum_file.close()
	# File with average conservation in 50 bins.
	with open(file_name + "_summary_bins.txt", "a+") as bins_file:
		if bins_file.read():
			bins_file.write("\n" + name)
		else:
			bins_file.write(name)
		for i in range(len(avg_phastcons)):
			if i < len(avg_phastcons):
				bins_file.write(str(avg_phastcons[i]) + "\t")
			else:
				bins_file.write(str(avg_phastcons[i]) + "\n")
	bins_file.close()
elif int(mode) == 2:
	avg_region = sum(region_data['phastcons']) / len(region_data['phastcons'])
	# Write the results in one file.
	# File with average conservation for the whole region.
	with open(file_name + "_summary_region.txt", "a+") as sum_file:
		if sum_file.read():
			sum_file.write("\n" + name + str(avg_region))
		else:
			sum_file.write(name + str(avg_region))
	sum_file.close()


