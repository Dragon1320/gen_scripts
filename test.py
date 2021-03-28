import copy

from ete3 import PhyloTree

from utils.orth import read_orthogroup_file
from utils.const import database_fp
from utils.misc import gene_to_sp
from utils.phylo import species_lookup, sp_id_to_name, root_tree, read_state_file, id_to_species, species_to_id, species_to_age, age_to_species, id_to_age, age_to_id, leaf_ids_to_names

focal_sp = "DMEL"
gene_threshold = 3

orthogroups = read_orthogroup_file(database_fp)

# NOTE: no need to filter out ogs that we use for the distance file
# those are guaranteed to be single copy until SLEB and were filtering out anything with <3 DMEL orthologs
potential_ogs = []

for og, genes in orthogroups.items():
  focal_genes = list(filter(lambda e: gene_to_sp(e) == focal_sp, genes))

  if len(focal_genes) >= gene_threshold:
    potential_ogs.append(og)

print("ogs: {}".format(len(potential_ogs)))

def gene_id_to_sp(tree_gene):
  return tree_gene[(tree_gene.rfind("_") + 1):]

for og in potential_ogs:
  genes = orthogroups[og]
  focal_genes = list(filter(lambda e: gene_to_sp(e) == focal_sp, genes))

  tree_fp = "./data/iqtree/group{}_alignment.fa.treefile".format(og)

  nw = None
  with open(tree_fp) as file:
    nw = file.readline().strip()

  tree = PhyloTree(nw, format = 1)

  # leaves = tree.get_leaves()

  # change all names to ids (above), easier to work with, as well as
  # ete3 needing 3 or less character long names for dating...
  # for leaf in leaves:
  #   leaf.name = species_to_id[leaf.name[(leaf.name.rfind("_") + 1):]]

  try:
    root_tree(tree)
  except:
    print("oof")

  # # cant always do this cos recursion stack overflow
  # try:
  #   tree = tree.collapse_lineage_specific_expansions()
  # except:
  #   # just for testing, actually analyse the tree when running this normally
  #   continue

  dmel_sp = "DMEL"
  dmel_id = species_to_id[dmel_sp]

  leaves = tree.get_leaves()
  times = set()

  for g1 in focal_genes:
    for g2 in focal_genes:

      ttree = copy.deepcopy(tree)
      l1 = ttree.search_nodes(name = g1.replace("|", "_"))[0]
      l2 = ttree.search_nodes(name = g2.replace("|", "_"))[0]

      tleaves = ttree.get_leaves()

      for leaf in tleaves:
        leaf.name = species_to_id[leaf.name[(leaf.name.rfind("_") + 1):]]

      event = ttree.get_common_ancestor(l1, l2)
      times.add(event.get_age(id_to_age))

  # for l1 in leaves:

  #   # i += 1
  #   # print("{}/{}".format(i, num_leaves))

  #   for l2 in leaves:
  #     if gene_id_to_sp(l1.name) != dmel_sp and gene_id_to_sp(l2.name) != dmel_sp:
  #       continue

  #     ttree = copy.deepcopy(tree)
  #     tl1 = ttree.search_nodes(name = l1.name)[0]
  #     tl2 = ttree.search_nodes(name = l2.name)[0]

  #     tleaves = ttree.get_leaves()

  #     for leaf in tleaves:
  #       leaf.name = species_to_id[leaf.name[(leaf.name.rfind("_") + 1):]]

  #     event = ttree.get_common_ancestor(tl1, tl2)
  #     times.add(event.get_age(id_to_age))

  print("*")

  rem = 0
  highest = 0

  for time in times:
    if time > 12:
      rem += 1

    if time > highest:
      highest = time

  # if highest >= 12:
  #   continue

  if len(times) - rem < 3:
    continue

  print(len(focal_genes), times, rem)
