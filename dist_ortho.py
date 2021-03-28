# given a list of fasta alignments, concat ones where all genes are in single copy and create a phylip file
from os import listdir
from os.path import isfile, join

from utils.fasta import read_fasta
from utils.misc import gene_to_sp

# const
num_species = 50

out_file_fp = "./out/infile"
alignment_dir = "./data/alignments"

seqs = {}

# onlyfiles lul
onlyfiles = [f for f in listdir(alignment_dir) if isfile(join(alignment_dir, f))]

for fname in onlyfiles:
  fp = join(alignment_dir, fname)

  fasta = read_fasta(fp)
  fasta_sp = [gene_to_sp(sp) for sp in fasta.keys()]

  # either duped or not all present
  if len(set(fasta_sp)) != len(fasta_sp) or len(fasta_sp) != 50:
    continue

  for gene, seq in fasta.items():
    sp = gene_to_sp(gene)

    if seqs.get(sp) == None:
      seqs[sp] = ""

    seqs[sp] += seq

# just phylip things - dont question it
outfile = open(out_file_fp, "w+")

outfile.write("{} {}\n".format(len(seqs), len(list(seqs.values())[0])))

for i, seq in seqs.items():
  outfile.write("{} {}\n".format(str(i).ljust(10), seq))

outfile.close()
