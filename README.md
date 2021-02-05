## data

All the required data for scripts to run.

- **/absense**
  - Bitscore, distances, and lengths files
- **/blast**
  - Self vs self blast for all species
- **/busco**
  - Dmel vs busco and busco v dmel blast files
- **/fasta**
  - Gene sequences in all species
- **drosophilaDatabase_neop**
  - Gene database (singletons, duplicable, rate)
- **orthogroups.txt**
  - Orthofinder orthogroups file
- **sequence_ids.txt**
  - Orthofinder gene ids file
- **species_ids.txt**
  - Orthofinder species ids file

## out

All the script outputs.

- **/absense**
  - Absense results
- **/gprofiler**
  - Gprofiler results
- **/plots**
  - Various plots from different scripts
- **busco_groups.txt**
  - Busco vs dmel reciprocal best hits -> orthogroups containing those genes
- **absense_example.png - temp**
- **paralogs.txt - temp**

## scripts

All the scripts.

- **busco**
  - Find all BUSCO Orthogroups
- **gp**
  - Runs GProfiler (fastevo, singleton, duplicable)
- **lengths**
  - Length comparison (fastevo, singleton, duplicable)
- **paralogs**
  - Paralog detection vs e-value cutoff (fastevo)
- **restriction**
  - Lineage restriction (fastevo)

## misc

- NOTE: not all scripts from the old repo have been refactored yet! - in particular the absesnse stuff
- TODO: ask about the boxplot stuff
- TODO: ask about the database
- TODO: list required libs somewhere (is there a python ver of package.json?)
