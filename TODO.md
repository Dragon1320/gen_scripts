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
- lengths:
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
