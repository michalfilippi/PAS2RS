#!/bin/bash

#
# Script for transformation of protein absolute structure to protein relative
# structure. Read README.md for more information.
#
# Usage: PAS2RS.sh input_file
# Usage: cat input_file | PAS2RS.sh
#

#create tmp file for coordinates
coordinates_file=$(tempfile)

# read stdin if no input file is given
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

# move coordinates to tmp file and pass it to the converter
cat "$input" > "$coordinates_file"
python3 protstr_converter.py "$coordinates_file"


#delete tmp file
rm "$coordinates_file"
