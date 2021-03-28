#!/bin/bash

function join_by { local IFS="$1"; shift; echo "$*"; }

# using these groups as theyre guaranteed to always be present
detected_group=("DMEL" "DSIM" "DSEC" "DERE" "DYAK" "DEUG" "DTAK" "DSUZ" "DBIA")

tmp_dir="./tmp"
nl=$'\n'

# absense script file
abs_fp="$1"
# distances file
dist_fp="$2"
# bitscores file
bit_fp="$3"
# gene lengths file
len_fp="$4"
# gene to run the analysis on (i dont wanna wait half a day for this to run)
gene="$5"

# prep the input files
rm -r "$tmp_dir"
mkdir -p "$tmp_dir/out"

include_only=`join_by , "${detected_group[@]}"`

bitscores_content="`cat "$bit_fp" | head -1`"$nl"`cat "$bit_fp" | grep "$gene"`"
lengths_content="`cat "$len_fp" | head -1`"$nl"`cat "$len_fp" | grep "$gene"`"
distances_content=`cat "$dist_fp"`

echo "$bitscores_content" > "$tmp_dir/bitscores.txt"
echo "$lengths_content" > "$tmp_dir/lengths.txt"
echo "$distances_content" > "$tmp_dir/distances.txt"

# run normally first
python3 "$abs_fp" --distfile "$tmp_dir/distances.txt" --scorefile "$tmp_dir/bitscores.txt" --genelenfile "$tmp_dir/lengths.txt" --predall True --out "$tmp_dir/out/normal"

# run ignoring everything other than the species specified
python3 "$abs_fp" --distfile "$tmp_dir/distances.txt" --scorefile "$tmp_dir/bitscores.txt" --genelenfile "$tmp_dir/lengths.txt" --includeonly "$include_only" --predall True --out "$tmp_dir/out/includeonly"

# calculate the spearman correlation
python3 ./abs_fit.py
