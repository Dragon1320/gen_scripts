import sys
from io import StringIO

from Bio import Phylo
from numpy import arange
from matplotlib import pyplot as plt

from utils.const import species_tree, database_fp, orthogroups_fp
from utils.data import get_fastevol_singletons
from utils.orth import find_gene_orthologs
from utils.misc import gene_to_sp

out_plt_fp = "./out/plots/plt_restriction.png"

tree = Phylo.read(StringIO(species_tree), "newick")
fsing = get_fastevol_singletons(database_fp)

data = []

for gene in fsing:
  orth = find_gene_orthologs(orthogroups_fp, gene)

  if len(orth) == 0:
    # TODO: ask what to do about this - theres orthogroups with only a single gene
    # so i wanna say these wont be genes without orthologs, but something that actually errored somewhere
    print("failed to find orthologs: ({})".format(gene))
    continue

  species = list(map(lambda e: gene_to_sp(e), orth))

  lowest_path = sys.maxsize
  lowest_mrca = None

  for sp in species:
    gene_sp = gene_to_sp(gene)

    mrca = tree.common_ancestor(sp, gene_sp)
    path = tree.get_path(mrca)

    lowest_path = len(path) if len(path) < lowest_path else lowest_path

    # if lowest_path was just set, set lowest_mrca
    if len(path) == lowest_path:
      lowest_mrca = mrca

  data.append(lowest_path)

# plot
plt.ylabel("num genes")
plt.xlabel("restriction")
plt.xticks(arange(12))

plt.hist(data, bins = arange(13), align = "left", edgecolor = "gray")

plt.savefig(out_plt_fp)

print(sorted(data))
