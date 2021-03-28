# given fasta files, generate a lengths file for given genes
from os import listdir
from os.path import isfile, join

from utils.misc import gene_to_sp
from utils.fasta import read_fasta
from utils.data import get_all_genes, get_fastevol_singletons
from utils.const import database_fp

# const
out_file_fp = "./out/lengths.txt"
fasta_dir = "./data/fasta"
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

outfile = open(out_file_fp, "w+")
outfile.write("GeneID\tGeneLength\n")

# get lengths from fasta files
for gene in final_genes:
  sp = gene_to_sp(gene)
  fasta_fp = "{}/{}.longest_only.faa".format(fasta_dir, sp)

  fasta = read_fasta(fasta_fp)
  seq = fasta[gene]

  outfile.write("{}\t{}\n".format(gene, len(seq)))

outfile.close()
