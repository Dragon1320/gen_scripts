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
