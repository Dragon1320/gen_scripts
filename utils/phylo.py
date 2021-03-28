from ete3.coretype.tree import TreeError

from utils.misc import gene_to_sp

pref_clades = [
  ["TCAS"],
  ["NVIT","CSOL","ACEP","MPHA","OBIR","BTER","MSAC","BTAB"],
  ["BMOR","TRNI","GMEL","HKAH","PPOL","PXUT","DPLE"],
  ["CNAS","AGAM","AALB"],
  ["SLEB"],
  ["DBUS","DGRI","DVIR","DNOV","DHYD","DNAV","DMOJ","DARI"],
  ["DWIL"],
  ["DOBS","DGUA","DPSE","DPER","DMIR"],
  ["DBIP","DANA"],
  ["DSER","DKIK"],
  ["DELE","DRHO"],
  ["DTAK","DSUZ","DBIA"],
  ["DEUG"],
  ["DERE","DYAK"],
]

def sp_id_to_name(node_name):
  return id_to_species[node_name]

def root_tree(tree):
  for clade in pref_clades:
    out_nodes = []

    for sp in clade:
      for node in tree.traverse():
        if sp in node.name:
          out_nodes.append(node)

    if len(out_nodes) == 0:
      continue

    if len(out_nodes) == 1:
      tree.set_outgroup(out_nodes[0])
      return
    else:
      root_node = tree.get_common_ancestor(out_nodes)

      try:
        tree.set_outgroup(root_node)
        return
      except TreeError:
        other_node = None

        for node in tree.get_leaves():
          if gene_to_sp(label_to_gene_name(node.name)) not in clade:
            other_node = node
            break
        else:
          continue

        tree.set_outgroup(other_node)
        root_node = tree.get_common_ancestor(out_nodes)
        tree.set_outgroup(root_node)
        return

  raise ValueError("could not find suitable root")

# utils - refactor/move later
def read_state_file(state_fp):
  nodes = {}

  with open(state_fp) as file:
    while True:
      line = file.readline().strip()

      # eof - idk how to actually handle eof in snek, tho this has always worked so lel
      if line == None or line == "":
        break

      # comments
      if line.startswith("#"):
        continue

      parts = line.split("\t")

      # header
      if parts[0] == "Node":
        continue

      node = parts[0]
      # site = parts[1]
      state = parts[2]

      if nodes.get(node) == None:
        nodes[node] = ""

      nodes[node] += state

  return nodes

species_lookup = [
  ("1", "DMEL", 1),
  ("2", "DSIM", 2),
  ("3", "DSEC", 2),
  ("4", "DERE", 3),
  ("5", "DYAK", 3),
  ("6", "DEUG", 4),
  ("7", "DTAK", 5),
  ("8", "DSUZ", 5),
  ("9", "DBIA", 5),
  ("10", "DELE", 6),
  ("11", "DRHO", 6),
  ("12", "DSER", 7),
  ("13", "DKIK", 7),
  ("14", "DBIP", 8),
  ("15", "DANA", 8),
  ("16", "DOBS", 9),
  ("17", "DGUA", 9),
  ("18", "DPSE", 9),
  ("19", "DPER", 9),
  ("20", "DMIR", 9),
  ("21", "DWIL", 10),
  ("22", "DGRI", 11),
  ("23", "DHYD", 11),
  ("24", "DNAV", 11),
  ("25", "DMOJ", 11),
  ("26", "DARI", 11),
  ("27", "DVIR", 11),
  ("28", "DNOV", 11),
  ("29", "DBUS", 11),
  ("30", "SLEB", 12),
  ("31", "CNAS", 13),
  ("32", "AGAM", 13),
  ("33", "AALB", 13),
  ("34", "BMOR", 14),
  ("35", "TRNI", 14),
  ("36", "GMEL", 14),
  ("37", "HKAH", 14),
  ("38", "PPOL", 14),
  ("39", "PXUT", 14),
  ("40", "DPLE", 14),
  ("41", "MSAC", 15),
  ("42", "BTAB", 15),
  ("43", "BTER", 15),
  ("44", "AMEL", 15),
  ("45", "OBIR", 15),
  ("46", "MPHA", 15),
  ("47", "ACEP", 15),
  ("48", "CSOL", 15),
  ("49", "NVIT", 15),
  ("50", "TCAS", 16),
]

id_to_species = dict(zip([i[0] for i in species_lookup], [i[1] for i in species_lookup]))
species_to_id = dict(zip([i[1] for i in species_lookup], [i[0] for i in species_lookup]))
species_to_age = dict(zip([i[1] for i in species_lookup], [i[2] for i in species_lookup]))
age_to_species = dict(zip([i[2] for i in species_lookup], [i[1] for i in species_lookup]))
id_to_age = dict(zip([i[0] for i in species_lookup], [i[2] for i in species_lookup]))
age_to_id = dict(zip([i[2] for i in species_lookup], [i[0] for i in species_lookup]))

def leaf_ids_to_names(tree):
  for leaf in tree.get_leaves():
    leaf.name = id_to_species[leaf.name]

def label_to_gene_name(label):
  return label.replace("_", "|")

def gene_name_to_label(gene):
  return gene.replace("|", "_")
