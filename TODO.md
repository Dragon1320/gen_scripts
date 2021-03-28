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

- *finish len graph
	- mostly done - needs title? (add after?), better labels
- *fix paralog script
	- there doesnt seem to be any issues with it?
	- is it that weird orthofinder stuff?
- *go through all denovo genes
	- some especially interesting ones:
		- LOC108005300
		- LOC108008871
		- LOC108016569
		- LOC108019732
- *make presence/absence trees nice
	- jesus that took forever
- better way to visualise gprofiler results
	- ensembl biomart? - think about this one a bit
- better way to visualise paralog results
	- id go for a heatmap like thing but we have 2 datasets...
- read through absense paper and make sure ik how to handle thingy
	- they do something like that, do we just replicate it with our data?
- make sure ik what im doing with paralog absense
	- were only looking at *that* part of the tree, right?
	- so how does one date the duplication events?

---

- (3) lengths file
- (2) look into (ete)? - start doing stuff for that - make sure i know what im doing here
- (2) visualisation of paralog shite
- (1) question in messages
- (1) start *something* on the writeup

- make /scripts reusable




- what does x mean in NodeA/x?

- so these either have paralogs or 0 genes for x species (what do?):
```
./data/iqtree/group924_alignment.fa
./data/iqtree/group6486_alignment.fa
./data/iqtree/group860_alignment.fa
./data/iqtree/group432_alignment.fa
./data/iqtree/group2302_alignment.fa
./data/iqtree/group201_alignment.fa
./data/iqtree/group751_alignment.fa
./data/iqtree/group901_alignment.fa
./data/iqtree/group6131_alignment.fa
./data/iqtree/group3731_alignment.fa
./data/iqtree/group1172_alignment.fa
```

- which ancestral nodes am i picking? all? only some?

- on the note of weird fucked up tree things - out of order and whatnot - were down to 107 seqs

- sets in python are weird (just like the string split escape bullshit)

---

- finish lengths file (or die trying)
- identify suitable orthogroups
- at ete3 example of how to analyse 1

---

what i did:
- lengths file just needs protdist ran
	- are we taking all (4) or only some (2) nodes?
	- theres very few nodes here...
- got ete3 stuff working
	- ???
- better paralog visualisation
	- am i doing it ok?

what i still need to do:
- distances file
- stuff in 'what i wanna get done by next week'
- that denovo gene stuff (lol)

what i wanna get done by next week:
- run absense
- that thing aoife mentioned - how good are the absense results actually (possibly?)

questions:
- questions in 'what i did'
- question in teams - still need it?
- not 100% sure what im doing here (maybe will be answered up until here?)
- what does x mean in NodeA/x?
- which ancestral nodes am i picking? all? only some?
- on the note of weird fucked up tree things - out of order and whatnot - were down to ~100 seqs

---

- TODO:
	- bitscores file
	- distances file
	- lengths file
	- heatmaps
	- how well absense predictions fit actual results

- NOTE:
	- ensure trees are rooted

- QUESTIONS:
	- so when rooting a tree, im guessing ete3 sometimes rearranges it in ways which add new nodes - cant get ancestral seq?
	- how well absense results fit actual hdf
	- presentation:
		- why its important
		- what has been done
		- what we did
		- what still needs to be done
		- do i include major things we tried that didnt work?
		- fitting everything into 10 minutes will be *interesting*

- colourblindness colours
- large fig text
- count for go terms in fast evolving singletons
- text overlap in absense fig
- duplicates singletons in paralog plots

- data
- gene dup slice
- bg project, fast evol sing came from?
- absense - how it works

---

- wrote the thing that tests absense predictions
	- am i doing it right?
- did heatmaps - work well (i think) just need some tweaks
- everything is (mostly) ready to run absense - just need the tree data and some small tweaks in my scripts
	- are we re-running ortholog absense with new genes for distances?
- currently working on the figures + presentation

- all figures that im putting in the writeup will basically be the same ones im including in the presentation

- figures:
	- larger text
	- better colours
	- absense ones:
		- nicer display
		- non-overlapping text
	- trees
		- version with full sp names
		- extend branches to right side?
		- highlight rate calculation/duplication status groups
		- a version with restriction levels
	- heatmaps
		- make sure the same colour in both means the same thing
		- colour-range display on the right
		- proper labels and all that stuff

presentation:
	- introduction
		- gene duplicability, its significance, etc.
		- basicaly what zoe had in her first few slides
	- data
		- the species tree
		- introduce the fast evolving singletons
		- orthofinder - quickly explain orthogroups (cos its a term ill probably be using a lot)
		- clearly state the aim of the project
	- legnths/go-terms
	- taxonomic restriction
	- absense (orthologs) - taxonomic restriction explained by hdf?
		- briefly mention de novo genes
	- absense (orthologs) w/ example results
	- absense (paralogs)
	- possibly emergency slides for questions im expecting to get (theres a very good chance ill get asked about de novo genes by aoife)
	- orthofinder doing weird stuff

tax restricted -> hdf idea -> absense to test it
only need to know methods

questions:
- into how much detail should i go in about the tools/methods - just mention stuff like orthofinder, 
- we need to re-run absense with new genes for distance calculations right? - if we want the 2 groups to be comparable which is what we intend
- should i mention that we tried BUSCO genes first and why that didnt work (briefly?)
- when will we have the data ready cos i ideally want to have the absense stuff ready to present on wednesday
- how should i mention the go terms in the presentation, would it be enough to say x group is enriched for y (or would a table on the slide in addition be better)
- i may need an explanation why orthofinder was being in case i want to mention it or im asked a question about it
- aim at start of slide 1 or somewhere in slide 2? - ooh i guess i can make just make it the title lel
- zoom link set to expire today

bannana genome venn diagram

---

- i was genuinely about to cry when i thought the distances were fucked up
- i checked out the bannana genome venn diagram LMAO (https://25.media.tumblr.com/tumblr_m723nmp5OM1qbh26io1_1280.jpg)
- absense has a paper out, only writing ~500 loc py
- i can def say im actually ok at python now
- cloning tree recursion error python being shit
- absense being annoying

colours:
- blue: #1f77b4
- orange: #ff7f0e
- green: #2ca02c (NOTE: dont mix with orange when distinction required)

TODO:
- absense no overlapping labels
- heatmaps
	- re-run with all data
	- make sure e-value labels are negatives

---

- more or less denovo genes in sing/fsing?
- how well absense results fit reality - how many/which ones should i run it on (if i want it to be correct)
- any better way to analyse the ortho/para absense results or is this the best we can do?
