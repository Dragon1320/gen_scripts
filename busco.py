from utils.blast import rec_best_hits
from utils.const import blast_busco_dmel_fp, blast_dmel_busco_fp, orthogroups_fp, database_fp
from utils.orth import find_gene_orthogroup

outfile_fp = "./out/busco_groups.txt"

hits = rec_best_hits(blast_dmel_busco_fp, blast_busco_dmel_fp)

hits_og = list(map(lambda e: find_gene_orthogroup(database_fp, e), hits))

file = open(outfile_fp, "w+") 

for og in hits_og:
  if og == None:
    continue

  file.write("{}\n".format(og))

file.close()
