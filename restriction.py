import sys
from io import StringIO

from Bio import Phylo
from numpy import arange, median
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
species_fp = "./data/species.txt"

species_abbr = {}

with open(species_fp) as f:
  lines = [l.strip() for l in f.readlines()]
  
  for l in lines:
    split = l.split("\t")

    species_abbr[split[0]] = split[1]

fsing = get_fastevol_singletons(database_fp)

print(len(fsing))

data = []

z = 0
zz = 0

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
  
  if lowest_path >= 4:
    z += 1

  if lowest_path != 0:
    zz += 1

  # debug
  lowest_mrca.color = "#1f77b4"

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
      sp.color = "#ff7f0e"
      undet_sp += 1

    undet_probs.append((sp.name, p))

  # debug
  if False:
  # if len(outside_sp) - undet_sp >= 2:
    leaves = tree.get_terminals()

    for leaf in leaves:
      dist_to_last_node = tree.distance(leaf) - leaf.branch_length
      leaf.branch_length = 0.75 - dist_to_last_node

    for sp, sp_name in species_abbr.items():  
      leaves = list(tree.find_elements(name = sp))

      if len(leaves) != 1:
        continue

      leaf = leaves[0]

      split = sp_name.split(" ")
      sp_abbr = "{}. {}".format(split[0][0], split[1])

      leaf.name = sp_abbr

    # thisll give an error the first time (on a machine without an display anyway) but thats fiiiiiine
    Phylo.draw(tree, branch_labels = lambda a: None)
    
    # plt.gca().invert_yaxis()
    fig = plt.gcf()
    fig.set_size_inches(10, 10)

    plt.ylabel("taxa", fontsize = 16, labelpad = 20)
    plt.xlabel("branch length", fontsize = 16, labelpad = 20)
    plt.tick_params(axis = "both", which = "major", labelsize = 12)

    plt.savefig("./out/trees/tree_test_{}_{}.png".format(gene[:gene.find("|")], len(outside_sp) - undet_sp))

  # print("{}: {} - {}".format(gene, lowest_path, list(map(lambda e: e[0], filter(lambda e: float(e[1]) >= 0.5, undet_probs)))))
  print("{}: {} - {}".format(gene, lowest_path, undet_probs))

print(z, zz)

# reset plt after all the tree stuff
plt.clf()

fig, ax = plt.subplots()

# plot
fig.set_size_inches(10, 10)

ax.set_xlabel("gene count", fontsize = 16, labelpad = 20)
ax.set_ylabel("restriction level", fontsize = 16, labelpad = 20)
ax.set_yticks(arange(12))
ax.tick_params(axis = "both", which = "major", labelsize = 12)

ax.hist(data, bins = arange(13), align = "left", edgecolor = "black", orientation = "horizontal")

line_handle = ax.axhline(median(data), color = "black", linestyle = "dashed", linewidth = 1)
ax.legend([line_handle], ["median"], fontsize = 12)

plt.savefig(out_plt_fp)

print(sorted(data))
