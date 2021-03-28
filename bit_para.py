# given gene trees, generate a lengths file

import copy
from functools import reduce

from ete3 import PhyloTree
from numpy import arange

from utils.orth import read_orthogroup_file
from utils.const import database_fp, gene_ids_fp, sp_ids_fp
from utils.misc import gene_to_sp, vk_lookup
from utils.phylo import species_lookup, sp_id_to_name, root_tree, read_state_file, id_to_species, species_to_id, species_to_age, age_to_species, id_to_age, age_to_id, leaf_ids_to_names, label_to_gene_name, gene_name_to_label
from utils.fasta import read_fasta
from utils.blast import read_blast_file

# const
focal_sp = "DMEL"
gene_threshold = 3
tree_dir = "./data/iqtree"
out_file_fp = "./out/bitscores_para.txt"
blast_dir = "./data/blast"

focal_sp_id = vk_lookup(sp_ids_fp, "{}.longest_only.faa".format(focal_sp))

orthogroups = read_orthogroup_file(database_fp)

# NOTE: no need to filter out ogs that we use for the distance file
# those are guaranteed to be single copy until SLEB and were filtering out anything with <3 DMEL orthologs
potential_ogs = []

for og, genes in orthogroups.items():
  focal_genes = list(filter(lambda e: gene_to_sp(e) == focal_sp, genes))

  if len(focal_genes) >= gene_threshold:
    potential_ogs.append(og)

print("ogs: {}".format(len(potential_ogs)))

final_genes = []

i = 0

for og in potential_ogs:

  # progress counter
  i += 1
  print("og: {} - {}/{}".format(og, i, len(potential_ogs)))

  genes = orthogroups[og]
  focal_genes = list(filter(lambda e: gene_to_sp(e) == focal_sp, genes))

  tree_fp = "{}/group{}_alignment.fa.treefile".format(tree_dir, og)

  nw = None
  with open(tree_fp) as file:
    nw = file.readline().strip()

  tree = PhyloTree(nw, format = 1)

  try:
    root_tree(tree)
  except:
    print("could not root tree: {}".format(og))
    continue

  for g1 in focal_genes:

    times = set()
    divergences = set()

    for g2 in focal_genes:
      # copy tree - we will need to edit the leaves to get divergence age
      # this sometimes gives a recursion error *sigh*
      try:
        ttree = copy.deepcopy(tree)

        l1 = ttree.search_nodes(name = gene_name_to_label(g1))[0]
        l2 = ttree.search_nodes(name = gene_name_to_label(g2))[0]

        # assign each leave an id based on divergence node
        tleaves = ttree.get_leaves()

        for leaf in tleaves:
          leaf.name = species_to_id[leaf.name[(leaf.name.rfind("_") + 1):]]

        event = ttree.get_common_ancestor(l1, l2)
        time = event.get_age(id_to_age)

        times.add(time)
        divergences.add((time, g2))
      except RecursionError:
        print("could not clone tree: {}".format(og))

    time_threshold = 3
    time_cutoff = 12
    above_cutoff = 0

    for time in times:
      if time > time_cutoff:
        above_cutoff += 1

    if len(times) - above_cutoff < time_threshold:
      continue

    found = False

    for gene, divs in final_genes:
      if gene == g1:
        found = True

        break

    if not found:
      final_genes.append((g1, divergences))

    print(g1, times, divergences)

sp_names = "\t".join(["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"])

outfile = open(out_file_fp, "w+")
outfile.write("None\t{}\n".format(sp_names))

i = 0

# get lengths from fasta files
for gene, divergences in final_genes:

  # progress counter
  i += 1
  print("gene {}/{}".format(i, len(final_genes)))

  gene_id = vk_lookup(gene_ids_fp, gene)

  bitscores = {}

  for n in arange(0, 12):
    bitscores[n] = []

  for t, subject in divergences:
    time = t - 1

    if time > 11:
      continue

    subject_id = vk_lookup(gene_ids_fp, subject)

    blast_fp = "{}/Blast{}_{}.txt.gz".format(blast_dir, focal_sp_id, focal_sp_id)
    blast = read_blast_file(blast_fp, gz = True)

    if bitscores.get(time) == None:
      bitscores[time] = []

    for entry in blast:
      query_blast = entry["query"]
      subject_blast = entry["subject"]

      if query_blast != gene_id or subject_blast != subject_id:
        continue

      bitscore = entry["score"]

      bitscores[time].append(bitscore)

  for time, scores in bitscores.items():
    if len(bitscores[time]) == 0:
      bitscores[time] = "0"
    else:
      bitscores[time] = str(reduce(lambda a, c: a + c, bitscores[time], 0) / len(bitscores[time]))

  outfile.write("{}\t{}\n".format(gene, "\t".join(bitscores.values())))

outfile.close()
