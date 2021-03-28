# given a list of blast files (focal sp vs everything else), generate a bitscore matrix for given genes
from os import listdir
from os.path import isfile, join
from functools import reduce

from utils.misc import gene_to_sp, vk_lookup
from utils.fasta import read_fasta
from utils.data import get_all_genes, get_fastevol_singletons
from utils.const import database_fp, species, sp_ids_fp, gene_ids_fp
from utils.blast import read_blast_file

# const
out_file_fp = "./out/bitscores.txt"
blast_dir = "./data/blast"
alignment_dir = "./data/alignments"

# genes = get_all_genes(database_fp)
genes = get_fastevol_singletons(database_fp)

# filter out ones used for distance calculation
dist_genes = []

onlyfiles = [f for f in listdir(alignment_dir) if isfile(join(alignment_dir, f))]

for fname in onlyfiles:
  fp = join(alignment_dir, fname)

  fasta = read_fasta(fp)
  fasta_sp = [gene_to_sp(sp) for sp in fasta.keys()]

  # either duped or not all present
  if len(set(fasta_sp)) != len(fasta_sp) or len(fasta_sp) != 50:
    continue

  dist_genes += fasta.keys()

final_genes = []

for db_gene in genes:
  if db_gene not in dist_genes:
    final_genes.append(db_gene)

sp_names = "\t".join(species)

outfile = open(out_file_fp, "w+")
outfile.write("None\t{}\n".format(sp_names))

i = 0

# get bitscores from blast files
for gene in final_genes:

  # progress counter
  i += 1
  print("gene {}/{}".format(i, len(final_genes)))

  gene_sp = gene_to_sp(gene)
  gene_id = vk_lookup(gene_ids_fp, gene)

  bitscores = {}

  for sp in species:
    sp_id = vk_lookup(sp_ids_fp, "{}.longest_only.faa".format(sp))
    gene_sp_id = vk_lookup(sp_ids_fp, "{}.longest_only.faa".format(gene_sp))

    blast_fp = "{}/Blast{}_{}.txt.gz".format(blast_dir, gene_sp_id, sp_id)
    blast = read_blast_file(blast_fp, gz = True)

    if bitscores.get(sp) == None:
      bitscores[sp] = []

    for entry in blast:
      query = entry["query"]

      if query != gene_id:
        continue

      bitscore = entry["score"]

      bitscores[sp].append(bitscore)

    if len(bitscores[sp]) == 0:
      bitscores[sp] = "0"
    else:
      bitscores[sp] = str(reduce(lambda a, c: a + c, bitscores[sp], 0) / len(bitscores[sp]))

  outfile.write("{}\t{}\n".format(gene, "\t".join(bitscores.values())))

outfile.close()
