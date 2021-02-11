import sys
from io import StringIO

from Bio import Phylo
from numpy import arange
from matplotlib import pyplot as plt

from utils.const import species_tree, database_fp, orthogroups_fp, database_fp, species_tree_fp
from utils.data import get_fastevol_singletons
from utils.orth import find_gene_orthologs
from utils.misc import gene_to_sp

# TODO: refactor
def undet_prob(gene, sp):
  species = []

  with open("./out/absense/Detection_failure_probabilities") as f:
    while True:
      l = f.readline().strip()

      if l.startswith("#"):
        continue

      if l == "" or l == None:
        break

      content = l.split("\t")

      # header line
      if content[0] == "Gene":
        species = content[1:]
        continue

      # non header line
      g = content[0]
      p = content[1:]

      if gene == g:
        idx = species.index(sp)

        return "0.0" if p[idx] == "Detected" else p[idx]

out_plt_fp = "./out/plots/plt_restriction.png"

fsing = get_fastevol_singletons(database_fp)

data = []

for gene in fsing:
  # instantiate the tree every time to remove any colouring
  tree = Phylo.read(species_tree_fp, "phyloxml")

  orth = find_gene_orthologs(database_fp, gene)

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

  # debug
  lowest_mrca.color = "cyan"

  # TODO: refactor
  # homology detection failure
  outside_sp = [e for e in tree.get_terminals() if e not in lowest_mrca.get_terminals()]

  undet_sp = 0
  undet_probs = []

  for sp in outside_sp:
    # a weird thing where if the most related ortholog is detected in TCAS,
    # the mrca will be the same as MSAC, etc.
    if lowest_path == 0:
      break

    p = undet_prob(gene, sp.name)

    if p == None:
      print("failed to find undetection probability: ({}, {})".format(gene, sp.name))
      continue

    if float(p) >= 0.5:
      sp.color = "red"
      undet_sp += 1

    undet_probs.append((sp.name, p))

  # debug
  if len(outside_sp) - undet_sp >= 2:
    # thisll give an error the first time (on a machine without an display anyway) but thats fiiiiiine
    Phylo.draw(tree, branch_labels = lambda a: None)
    
    # plt.gca().invert_yaxis()
    fig = plt.gcf()
    fig.set_size_inches(10, 10)

    plt.savefig("./out/trees/tree_test_{}.png".format(gene[:gene.find("|")]))

  # print("{}: {} - {}".format(gene, lowest_path, list(map(lambda e: e[0], filter(lambda e: float(e[1]) >= 0.5, undet_probs)))))
  print("{}: {} - {}".format(gene, lowest_path, undet_probs))

# reset plt after all the tree stuff
plt.clf()

# plot
plt.ylabel("gene count")
plt.xlabel("restriction level")
plt.xticks(arange(12))

plt.hist(data, bins = arange(13), align = "left", edgecolor = "gray")

plt.savefig(out_plt_fp)

print(sorted(data))
