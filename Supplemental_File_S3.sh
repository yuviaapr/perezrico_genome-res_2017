#!/bin/bash

# Script to generate files with the proportion of specific regions (e.g. gene bodies) overlapping with typical enhancers, super-enhancers and bootstrap controls.
# The first line in the output files corresponds to typical enhancers or super-enhancers and the following n lines to n bootstrap control iterations.
# This script relies on bedtools shuffle and bedtools annotate.
# Usage: ./Supplemental_File_S3.sh regions.bed typical_enhancers.bed super_enhancers.bed iterations_bootstrap chrom_sizes.txt prefix
# 	regions.bed = BED file containing the coordinates of genomic features.
#	typical_enhancers.bed = BED file containing the coordinates of typical enhancers.
#	super_enhancers.bed = BED file containing the coordinates of super-enhancers.
# 	iterations_bootstrap = number of iterations for bootstrap analyses.
#	chrom_sizes = file with chromosome sizes.
#	prefix = prefix for the output files.

#Initialize variables.

regions_f=$1
typical=$2
super_enhancer=$3
iterations_n=$4
chr_sizes=$5
pre_name=$6
n=1

# Assign names to the output files.

typical_file=${pre_name}_$( echo -n typical_enhancers.txt )
super_file=${pre_name}_$( echo -n super_enhancers.txt )

# Calculate proportions for typical enhancers and super-enhancers.

bedtools annotate -counts -i $regions_f -files $typical | awk 'BEGIN{cnt=0};{if($4!=0){cnt=cnt+1}};END{print"Typical_enhancers\t"cnt/NR}' > $typical_file

bedtools annotate -counts -i $regions_f -files $super_enhancer | awk 'BEGIN{cnt=0};{if($4!=0){cnt=cnt+1}};END{print"Super_enhancers\t"cnt/NR}' > $super_file

# Iterate the analysis using shuffled control regions.

while (( $n <= $iterations_n )); do
	# Typical enhancers
	bedtools shuffle -i $typical -g $chr_sizes -chrom -noOverlapping | bedtools annotate -counts -i $regions_f -files stdin | awk 'BEGIN{cnt=0};{if($4!=0){cnt=cnt+1}};END{print"control\t"cnt/NR}' >> $typical_file
	# Super-enhancers
	bedtools shuffle -i $super_enhancer -g $chr_sizes -chrom -noOverlapping | bedtools annotate -counts -i $regions_f -files stdin | awk 'BEGIN{cnt=0};{if($4!=0){cnt=cnt+1}};END{print"control\t"cnt/NR}' >> $super_file
	n=$(( n+1 ))	
done


