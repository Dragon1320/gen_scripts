compare lengths - mann whitney u test
fast evolving signletons Vs everything else
two tailed test

assign how taxonomically restricted

suzikii upwards
lengths from fasta file + take average
compare lengths of fast evolving singletons vs everything else

---

file with bitscores and with distances
busco genes - conserved, generally single copy
DSUZ as focal first, then maybe DMEL focal species, better for finding busco genes

species ids
sequence ids
find orthologs in all other species (orthogroups) then find bitscores

if multiple orthologs, take avg bitscore

length file optional

---

blast vs genome
reciprocal best hits
0.001
melanogaster
orthogroups
align
convert to phylyp format

- prediction intervals

---

q:
- missing ids (verify)
- BUSCO: ancestral vs ancestral_variants

w:
- blastp BUSCO genes vs DMEL
- take all unique DMEL gene hits where e < 0.001
- find orthologs where at least 1 gene exists in all other drosophila species
- align, phylyp, distances

---

TODO:
- ~~gen final dist file~~
- ~~fix and re-run lengths/bitscore (gene amount needs to match)~~

SANITY_CHECKS:
- for lengths and bitscores - using all genes in db where S/D and excluded reason is null

QUESTIONS:
- im asuming this is fine (blast): CFastaReader: Hyphens are invalid and will be ignored around line 12318, 14176, 26192 (https://www.biostars.org/p/394195/)
- some fasta files have \*, muscle doesnt seem to like those - stop codons, what do? (https://www.biostars.org/p/148084/)
- bitscore.py - there are some genes in orthogroups.txt that dont have corresponding ids in sequenceids.txt
	- Dyak\EG_EG0007.10|XP_002100128.1|DYAK (orth)
	- Dyak\EG:EG0007.10|XP_002100128.1|DYAK (idfile)
	- l_1_G0289|NP_572621.1|DMEL (orth)
	- l(1)G0289|NP_572621.1|DMEL (idfile)
- database lengths? - do we need them? what are they?
- how to do graphs/data vis - recommendations?
- phd

RAMBLINGS:
- the seemingly endless pain with \/ and \
- things just like to replace things for no apparent reason

PROCEDURE:
- lengths:
	- get lengths from fasta file
- bitscores:
	- get orthologs for each gene
	- check all blast files for gene vs ortholog bitscore
	- if multiple orthologs for gene in x species, take avg
- dists:
	- reciprocal best hits eukaryota vs DMEL
	- in those DMEL genes, reciprocal best hits vs all other sp
	- take only genes where there is an ortholog found in all sp (by method above)
	- make alignments for all these orthologs per busco gene
	- concat genes from alignments (but this time by species) and put it phylip format
	- protdist default settings, non-interleaved

GOALS:
- fix any issues
- run absense
- analyse le data

TODO:
- how to analyse the data
-	adapt for paralogs

---

QUESTIONS:
- is it fine lol
- stuff from last week?
- are we trying other focal species?
	- can easily clean up my scripts to automate that
	- will need to essentially redo what was done for the db tho
- a bit on displaying the data

THOUGHTS:
- interesting res.png - homolog detection failure possible explanation - how likely? - look at what the paper was trying to prove again!
- i guess if we used DMEL as the focal species, wed be able to examine the really tr genes
- the absense scripts rly dont look too bad, only a couple 100 loc
- i def need to have more of a think about the data and look through it

---

NEW:
- clean up repo + document stuff
- detection failure vs lineage restriction

TODO:
- blast agaionst itserlf for fast ev singletons
- characterise gene functions for singletons
- functional profile of singletons - similar to duplicable genes?
	- gprofiler - need to play around with gene ids - might need flybase ids?
	- melanogaster orghologs
	- muiltiple test6iong correction: false discovery rate (less strict than bonferroni - lil strict for func enrichment)
	- fastr evol, other ones, duplicable ones - background: all genes (in genome)
	- significant depletions (opt?)









---

ref: https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/

- is what i did for trg correct? stats? - need more specifics!
- eval cutoff - saw e = 0.1 in papers (https://academic.oup.com/mbe/article/35/1/107/4554431#113623311)
- 'CFastaReader: Hyphens are invalid and will be ignored around line 12318' error again



---

latest in tree where not detected, see if its detected in abs
incrasing e val for ingroup


iqtree - ancestral seq - busco groups
- >=2 independent dups
- send results

---

- need an explanation of the qtree stuff cos im dumb
- how did we end up with the database stuff
- am i doing the paralog stuff correctly, cutoffs to use
- the wednesday talks

---

lengths other dup AND sing
enrichment (better groups - use database)
take non stat sig ones as well for fastevo
	- look how many geens have these terms assigned

histogram of restriction levels
lineage specific gene loss

group fastevol
- restricted or not
	- restricted - denovo vs hdf

---

- finish len graph
	- mostly done - needs title? (add after?), better labels
- fix paralog script
	- there doesnt seem to be any issues with it?
	- is it that weird orthofinder stuff?
- go through all denovo genes
	- some especially interesting ones:
		- LOC108005300
		- LOC108008871
		- LOC108016569
		- LOC108019732
- make presence/absence trees nice
- read through absense paper and make sure ik how to handle thingy
- better way to visualise gprofiler results
	- ensembl biomart?
- better way to visualise paralog results
	- id go for a heatmap like thing but we have 2 datasets...
- make sure ik what im doing with paralog absense





















