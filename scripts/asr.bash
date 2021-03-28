#! /bin/bash

while IFS= read -r line; do
	alignment="group"$line"_alignment.fa"

	echo "$alignment"
	iqtree -s "$alignment" -bb 1000 -nt 2 -mset WAG,LG,JTT -asr
done < "$1"
