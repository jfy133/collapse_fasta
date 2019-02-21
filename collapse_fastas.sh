#!/usr/bin/env bash

## Author: James A. Fellows Yates (jfy133[at]gmail.com)
## Date: February 2019

usage="
	Usage: 
	collapse_fasta.sh <output_name> /<path_to/<output_dir>/ /<path_to>/<fasta_1>.fa /<path_to>/<fasta_2>.fa [...]
		
	Description:
		Takes as input: a name, an output directory and a list of FASTA files. 
		It combines all separate files into a single file with a single header 
		(specified by the name argument), and makes a-tab separated coordinates 
		file with start/end coordinates of each FASTA entry (i.e. reference 
		genome, chromosomes and/or contigs) in the new collapsed fasta.

	Notes: 
	  * FASTAs can be a mixture of uncompressed or gzipped (.gz).
	  * Output columns are: 
	    1) input_file 
	    2) fasta_entry_header 
	    3) fasta_entry_length 
	    4) start_coordinate 
	    5) end_coordinate

	Options:
		--help / -h 	show this help text
		"

## Error messages
if [[ "$1" == "-h" ]]; then
	echo "$usage"
	exit 0
elif [[ "$1" == "--help" ]]; then
	echo "$usage"
	exit 0
elif [[ $# -le 1 ]]; then
	echo "
	Error: Incorrect input arguments! 

	$usage"
	exit 1
elif [[ -d $1 ]]; then
	echo "
	Error: Incorrect input arguments! 

	$usage"
	exit 1
elif [[ -f $1 ]]; then
	echo "
	Error: Incorrect input arguments! 

	$usage"
	exit 1
elif [[ -f $2 ]]; then
	echo "
	Error: Incorrect input arguments! 

	$usage"
	exit 1
elif [[ ! -d $2 ]]; then
	echo "
	Error: Incorrect input arguments! 

	$usage"
	exit 1
fi

## Set variables, slicing array to remove out_dir from file list and find if 
## invalid extension mix
name=$1
out_dir="$(readlink -f "$2")"
file_array=("${@:1}")
file_array=("${file_array[@]:2}")

## Removing old files
#echo "...[If run previously] removing old temporary files..."

if [[ -f "$out_dir"/temp_collapsed_"$name".coords ]]; then
	rm "$out_dir"/temp_collapsed_"$name".coords
fi

if [[ -f "$out_dir"/temp_collapsed_"$name".fa ]]; then
	rm "$out_dir"/temp_collapsed_"$name".fa
fi

if [[ -f "$out_dir"/collapsed_"$name".fa ]]; then
	echo "
	Error: $out_dir/collapsed_$name.fa already exists!

	Please remove the old file or provide a new output name.
	"
	exit 1
fi

if [[ -f "$out_dir"/collapsed_"$name".coords ]]; then
	echo "
	Error: $out_dir/collapsed_$name.coords already exists!

	Please remove the old file or provide a new output name.
	"
	exit 1
fi



## Generate per-fasta coordinate system
## Combine fastas into a single file (currently with headers)
## awk oneliner: https://www.danielecook.com/generate-fasta-sequence-lengths/

#echo "...generating contig coordinates and generating collapsed FASTA..."

for i in "${file_array[@]}"; do
	file_name="$i"
	if [[ "${i##*.}" == 'gz' ]]; then
		awk -v file_name="$file_name" \
		'$0 ~ /^>/ {
			print c; c=0;printf file_name "\t" substr($0,2) "\t"; 
		} $0 !~ ">" {
			c+=length($0);
		} END {
			print c; 
		}' <(zcat "$i") >> "$out_dir"/temp_collapsed_"$name".coords
		zcat "$i" >> "$out_dir"/temp_collapsed_"$name".fa
	else
		awk -v file_name="$file_name" \
		'$0 ~ /^>/ {
			print c; c=0;printf file_name "\t" substr($0,2) "\t"; 
		} $0 !~ ">" {c+=length($0);} END {
			print c; 
		}' "$i" >> "$out_dir"/temp_collapsed_"$name".coords
		cat "$i" >> "$out_dir"/temp_collapsed_"$name".fa
	fi
done

## Clean empty lines
sed -i '/^$/d' "$out_dir"/temp_collapsed_"$name".coords

if [[ ! -f "$out_dir"/temp_collapsed_"$name".coords ]]; then
	echo "Oops... something went wrong..."
	exit 1
fi

if [[ ! -f "$out_dir"/temp_collapsed_"$name".fa ]]; then
	echo "Oops... something went wrong..."
	exit 1
fi

#echo "...calculating new FASTA coordinates..."

awk -F '\t' '{sum += $3; start=(sum-$3)+1} {print $0"\t"start"\t"sum}' \
"$out_dir"/temp_collapsed_"$name".coords > "$out_dir"/collapsed_"$name".coords

#echo "...stripping FASTA headers..."

## From: https://www.unix.com/302533338-post2.html
awk '/>/&&c++>0 {
	next
} 1 ' "$out_dir"/temp_collapsed_"$name".fa | awk '!/^>/ { 
	printf "%s", $0; n = "\n" 
} /^>/ { 
	print n $0; n = "" 
} END { 
	printf "%s", n 
}' | fold -w 80 > "$out_dir"/collapsed_"$name".fa
sed -i "s/>.*/>$name/g" "$out_dir"/collapsed_"$name".fa

#echo "...cleaning up..."

rm "$out_dir"/temp_collapsed_"$name".coords
rm "$out_dir"/temp_collapsed_"$name".fa

#echo "...Done!"

exit 0