#! /bin/bash

while IFS= read -r line; do
	alignment="group"$line"_alignment"

	find . -name "$alignment*.*" -exec cp '{}' ../busco/ \;
done < "$1"
