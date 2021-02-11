from utils.blast import read_blast_file, best_hits, rec_best_hits
from utils.const import blast_busco_dmel_fp, blast_dmel_busco_fp, orthogroups_fp, database_fp
from utils.orth import find_gene_orthologs, read_orthogroup_file

res = read_orthogroup_file(database_fp)

print(res)
