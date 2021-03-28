from Bio import Phylo
from matplotlib import pyplot as plt

tree_fp = "./data/tree.xml"
species_fp = "./data/species.txt"

species = {}

with open(species_fp) as f:
  lines = [l.strip() for l in f.readlines()]
  
  for l in lines:
    split = l.split("\t")

    species[split[0]] = split[1]

tree = Phylo.read(tree_fp, "phyloxml")

leaves = tree.get_terminals()

rate_nodes = []
dup_nodes = []

for leaf in leaves:
  dist_to_last_node = tree.distance(leaf) - leaf.branch_length
  leaf.branch_length = 0.75 - dist_to_last_node

  if leaf.name in ["DEUG", "DTAK", "DSUZ", "DBIA"]:
    rate_nodes.append(leaf)

  if leaf.name in ["DMEL", "DSIM", "DSEC", "DERE", "DYAK"]:
    dup_nodes.append(leaf)

rate_anc = tree.common_ancestor(*rate_nodes)
dup_anc = tree.common_ancestor(*dup_nodes)

rate_anc.color = "#ff7f0e"
dup_anc.color = "#1f77b4"

for sp, sp_name in species.items():  
  leaves = list(tree.find_elements(name = sp))

  if len(leaves) != 1:
    continue

  leaf = leaves[0]

  split = sp_name.split(" ")
  sp_abbr = "{}. {}".format(split[0][0], split[1])

  leaf.name = sp_abbr

# plot
Phylo.draw(tree, branch_labels = lambda a: None)

fig = plt.gcf()
fig.set_size_inches(10, 10)

plt.ylabel("taxa", fontsize = 16, labelpad = 20)
plt.xlabel("branch length", fontsize = 16, labelpad = 20)
plt.tick_params(axis = "both", which = "major", labelsize = 12)

plt.savefig("./out/trees/species_highlight.png")
