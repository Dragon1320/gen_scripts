from os import listdir
from os.path import isfile, join

from Bio import Phylo
from matplotlib import pyplot as plt
from ete3 import PhyloTree
from ete3.coretype.tree import TreeError

from utils.const import ingroup_sp, blast_busco_dmel_fp, blast_dmel_busco_fp, database_fp, translation_fasta_dir
from utils.orth import find_gene_orthogroup
from utils.blast import rec_best_hits
from utils.fasta import read_fasta
from utils.phylo import species_lookup, sp_id_to_name, root_tree, read_state_file, id_to_species, species_to_id, species_to_age, age_to_species, id_to_age, age_to_id, leaf_ids_to_names

# const
out_file_fp = "./out/infile"
tree_dir = "./data/iqtree3"

# onlyfiles lul
onlyfiles = [f for f in listdir(tree_dir) if isfile(join(tree_dir, f)) and f.endswith(".state")]

### DMEL orthologs
# og_hits = {}

# hits = rec_best_hits(blast_dmel_busco_fp, blast_busco_dmel_fp)

# for hit in hits:
#   og_hits[find_gene_orthogroup(database_fp, hit)] = hit

### get busco groups
groups = []

for f in onlyfiles:
  og = f[:(f.find("."))]

  groups.append(og)

# with open("./data/groups2.txt") as f:
#   groups = list(map(lambda e: e.strip(), f.readlines()))

seqs = {}

c = 0
cc = 0

### if, by some chance, there is a god (or gods - i dont judge), im definitely going to some special kind of hell for writing this code
for group in groups:

  # progress counter
  c += 1
  print("group: {} - {}/{}".format(group, c, len(groups)))

  # if group != "OG0005937":
  #   continue

  # temp - clean up later
  # data_path = "./data/iqtree3"
  tree_fp = "{}/{}.out.fa.treefile".format(tree_dir, group)

  nw = None
  with open(tree_fp) as file:
    nw = file.readline().strip()

  tree = PhyloTree(nw, format = 1)

  leaves = tree.get_leaves()

  try:
    root_tree(tree)
  except:
    print("could not root tree: {}".format(group))
    continue

  # change all names to ids (above), easier to work with, as well as
  # ete3 needing 3 or less character long names for dating...
  for leaf in leaves:
    leaf.name = species_to_id[leaf.name[(leaf.name.rfind("_") + 1):]]

  tree = tree.collapse_lineage_specific_expansions()

  # split nodes based on which ancestral node they belong to
  sp_groups = []
  last_age = 0

  for i, sp, age in species_lookup:
    if age > last_age:
      last_age = age
      sp_groups.append([(i, sp, age)])

      continue

    sp_groups[-1].append((i, sp, age))

  # make sure that the ancestral node of DMEL and all nodes in a group is constant
  # TODO: make sure theres only one of these!
  dmel_nodes = tree.search_nodes(name = species_to_id["DMEL"])
  if len(dmel_nodes) != 1:
    continue

  dmel_node = dmel_nodes[0]

  ancestral_nodes = []

  # measures how many clades which should give the same ancestral node dont
  out_of_order = 0
  # 2 is quite good, 1 for testing shite
  out_of_order_allowed = 666

  for species in sp_groups:
    nodes = []

    for i, sp, age in species:
      sp_nodes = tree.search_nodes(name = i)

      nodes += sp_nodes

    ancestors = set()

    for node in nodes:
      ancestor = tree.get_common_ancestor(dmel_node, node)

      ancestors.add(ancestor)

    if len(ancestors) == 0:
      break

    if len(ancestors) != 1:

      out_of_order += 1

      ancestors = { tree.get_common_ancestor(*ancestors) }

      # break

    # print(ancestors, species)

    ancestral_nodes.append(list(ancestors)[0])

  if out_of_order > out_of_order_allowed:
    continue

  # make sure the ancestral node ages are actually in order
  ordered = True

  # exclude DMEL cos that can be a bit weird in the tree sometimes cos reasons
  # if everything else is in order, it should be fiiiiine
  for i, node in enumerate(ancestral_nodes):
    # no next node
    if i == len(ancestral_nodes) - 1:
      break

    a = node
    b = ancestral_nodes[i + 1]

    # b should be the older node so if we get the common ancestor of both,
    # it should always be b
    anc = tree.get_common_ancestor(a, b)
    # print(anc.name, b.name)
    # print(anc.get_leaves())

    if anc != b:
      ordered = False
      break

  if not ordered:
    continue

  node_count = 12

  # only take the first node_count nodes
  ancestral_nodes = ancestral_nodes[:node_count]

  # node_count - number of ancestral nodes which need to be in order
  # NOTE: if there isnt enough, data, set this to 10 instead
  if len(ancestral_nodes) < node_count:
    continue

  # TODO: check if we can leave this out - we lose a lot of data by doing this
  # filter out groups if ancestral nodes are not unique
  if len(ancestral_nodes) != len(set(ancestral_nodes)):
    # leaf_ids_to_names(tree)
    # print(tree)
    # exit(1)

    # temp
    # if len(ancestral_nodes) - len(set(ancestral_nodes)) <= 1:
    #   c += 1

    #   print(len(ancestral_nodes), len(set(ancestral_nodes)))
    #   print(tree)

    #   continue

    # print(len(ancestral_nodes), len(set(ancestral_nodes)))

    continue

  # c += 1
  # print(c)

  # continue

  # make sure all nodes have an ancestral node labelled assigned
  # my guess is some disappear after rooting the tree
  found = True
  
  for i, node in enumerate(ancestral_nodes[1:node_count]):
    if node.name == "":
      found = False
      break

  if not found:
    continue

  # handle dmel genes
  dmel_phylip_id = 0

  fasta_file_fp = "{}/{}.out.fa".format(tree_dir, group)
  fasta_file = read_fasta(fasta_file_fp)

  # dmel_id = og_hits.get(int(group))
  # if dmel_id == None:
  #   print("wut")
  #   exit(1)

  # all genes are guaranteed to be singletons (up until SLEB at least)
  dmel_ids = [gene for gene in fasta_file.keys() if gene[-4:] == "DMEL"]

  # just in case
  if len(dmel_ids) != 1:
    print("wut")
    exit(1)

  dmel_seq = fasta_file.get(dmel_ids[0])
  if dmel_seq == None:
    print("oof")
    exit(1)

  if seqs.get(dmel_phylip_id) == None:
    seqs[dmel_phylip_id] = ""

  seqs[dmel_phylip_id] += dmel_seq

  # minus DMEL - handled separately
  for i, node in enumerate(ancestral_nodes[1:node_count]):
    slash_idx = node.name.find("/")
    ancestral_node_name = node.name[:(None if slash_idx == -1 else slash_idx)]
    # cos we skipped dmel, which will be 0
    idx = i + 1

    state_file_fp = "{}/{}.out.fa.state".format(tree_dir, group)
    state = read_state_file(state_file_fp)

    if seqs.get(idx) == None:
      seqs[idx] = ""

    seqs[idx] += state[ancestral_node_name]

  cc += 1

print(cc)
# exit(1)

# just phylip things - dont question it
outfileFp = "./out/infile"
outfile = open(outfileFp, "w+")

outfile.write("{} {}\n".format(len(seqs), len(list(seqs.values())[0])))

for i, seq in seqs.items():
  outfile.write("{} {}\n".format(str(i).ljust(10), seq))

outfile.close()
