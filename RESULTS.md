# results

- link: 

## a

- all genes taken from database
- two-sided mannwhitneyu test - cmp length of fast evolving singletons (where rate > mean + std * 2) vs everything else
- MannwhitneyuResult(statistic=53696.5, pvalue=3.062168921976593e-08)

## b

- gprofiler results
- file names: gprofiler_<group>_<mode>
  - group:
    - fastevol - fast evolving singletons
    - single - all other singletons
    - other - non-singletons
  - mode:
    - normal - measure_underrepresentation = False
    - underr - measure_underrepresentation = True

## c

- absense results

## d

- taxonomic restriction and detection probability in fast evolving singletons
- format:
  - gene - gene name
  - restriction - how taxonomically restricted the gene is
  - undetected - number of species in which ortholog detection failure probability >= 0.5
  - detected - number of species in which ortholog detection failure probability < 0.5
