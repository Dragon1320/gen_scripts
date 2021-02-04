import sys
from io import StringIO
from functools import reduce
from enum import Enum

from Bio import Phylo

from utils.db_init import get_db_cursor

# const
orthogroupPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"

def undet_prob(gene, sp):
  species = []

  with open("./abs_res/Detection_failure_probabilities") as f:
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

# const
kDbPath = "../drosophilaDatabase_neop"

species_tree = "(TCAS,(((((BMOR,TRNI),(GMEL,HKAH)),((PPOL,PXUT),DPLE)),((CNAS,(AGAM,AALB)),(SLEB,(((DGRI,((DHYD,(DNAV,(DMOJ,DARI))),(DVIR,DNOV))),DBUS),(DWIL,(((DOBS,DGUA),((DPSE,DPER),DMIR)),((DBIP,DANA),((DSER,DKIK),((DELE,DRHO),((DTAK,(DSUZ,DBIA)),(DEUG,((DERE,DYAK),(DMEL,(DSIM,DSEC)))))))))))))),((MSAC,BTAB),(((BTER,AMEL),(OBIR,(MPHA,ACEP))),(CSOL,NVIT)))));"

# setup
handle = StringIO(species_tree)
tree = Phylo.read(handle, "newick")

# c = get_db_cursor(kDbPath)

# temp
# Phylo.draw_ascii(tree)

# db_temp = ["CG33687|NP_001027130.2|DMEL"]

# c.execute("SELECT id FROM processed_trees")
# db_genes = list(map(lambda gene: gene["id"], c.fetchall()))

# shork

# db_genes = []

# with open("./abs_res/Detection_failure_probabilities") as f:
#   while True:
#     l = f.readline().strip()

#     if l.startswith("#"):
#       continue

#     if l is "":
#       break

#     sp = l.split("\t")

#     if sp[0].find("DSUZ") == -1:
#       continue
    
#     db_genes.append(sp[0])

# a

# 1 gene can only be in 1 orthogroup (load ALL genes in file and compare len of list vs set)

# print("total genes: {}".format(len(db_genes)))

# curr = 0

def find_orthogroup_genes(path, gene):
  with open(path) as f:
    while True:
      l = f.readline()

      if l == "" or l == None:
        break

      genes = l[(l.find(":") + 1):].strip().split(" ")

      if gene in genes:
        return genes

  return None

cur = 0

from utils.genes import get_fastevol, get_single

db_genes = get_fastevol()

print("total genes: {}".format(len(db_genes)))

# tempa = []
# tempb = []

# temp - testing
# db_genes = ["LOC108009337|XP_016929104.1|DSUZ"]

out_file = open("./trg_abs_res.txt", "w+")
out_file.write("{}\t{}\t{}\t{}\n".format("gene", "restriction", "undetected", "detected"))

for gene in db_genes:
  orth = find_orthogroup_genes(orthogroupPath, gene)

  if orth == None:
    print("failed to find orthologs: ({})".format(gene))
    continue

  species = list(map(lambda e: e[(e.rfind("|") + 1):], orth))

  lowest_path = sys.maxsize
  lowest_mrca = None

  for sp in species:
    gene_sp = gene[(gene.rfind("|") + 1):]

    mrca = tree.common_ancestor(sp, gene_sp)
    path = tree.get_path(mrca)

    lowest_path = len(path) if len(path) < lowest_path else lowest_path

    # if lowest_path was just set, set lowest_mrca
    if len(path) == lowest_path:
      lowest_mrca = mrca

  outside_sp = list(map(lambda e: e.name, [e for e in tree.get_terminals() if e not in lowest_mrca.get_terminals()]))
  # outside_orth = list(map(lambda sp: next(e for e in orth if e[(e.rfind("|") + 1):]), outside_sp))

  undet_sp = 0
  undet_probs = []

  for sp in outside_sp:
    # a weird thing where if the most related ortholog is detected in TCAS,
    # the mrca will be the same as MSAC, etc.
    if lowest_path == 0:
      break

    p = undet_prob(gene, sp)

    if p == None:
      print("failed to find undetection probability: ({}, {})".format(gene, sp))
      continue

    if float(p) >= 0.5:
      undet_sp += 1

    undet_probs.append((sp, p))

  print("{}. {}: {} - {}".format(cur, gene, lowest_path, list(map(lambda e: e[1], undet_probs))))

  out_file.write("{}\t{}\t{}\t{}\n".format(gene, lowest_path, undet_sp, len(outside_sp) - undet_sp))
  # out_file.write("{}: {} - {}".format(gene, lowest_path, list(map(lambda e: e[2], undet_probs))))
  # out_file.write("{}: {}".format(gene, lowest_path))

  cur += 1

out_file.close()

# from matplotlib import pyplot
# from scipy.stats import pearsonr, spearmanr

# pyplot.scatter(tempa, tempb)
# pyplot.savefig("trg_abs_plot_test.png")

# pcor, _ = pearsonr(tempa, tempb)
# scor, _ = spearmanr(tempa, tempb)

# print("pearson: {}\nspearman: {}".format(pcor, scor))
