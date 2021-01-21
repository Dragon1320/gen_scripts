import sys
from io import StringIO
from functools import reduce
from enum import Enum

from Bio import Phylo

from utils.db_init import get_db_cursor

# True - homologs of the gene would be detected in all species
# False - homologs of the gene may not be detected in at least 1 species
def homolog_det(gene):
  with open("./abs_res/Detection_failure_probabilities") as f:
    while True:
      l = f.readline().strip()

      if l.startswith("#"):
        continue

      if l is "":
        break

      sp = l.split("\t")

      if sp[0] != gene:
        continue

      for ip in range(1, len(sp)):
        if sp[ip] == "Detected" or sp[ip] == "0.0":
          continue

        return False

      return True

  return False

def homolog_det_count(gene):
  with open("./abs_res/Detection_failure_probabilities") as f:
    detected = 0

    while True:
      l = f.readline().strip()

      if l.startswith("#"):
        continue

      if l is "":
        break

      sp = l.split("\t")

      if sp[0] != gene:
        continue

      for ip in range(1, len(sp)):
        if sp[ip] == "Detected" or sp[ip] == "0.0":
          detected += 1

      return detected

  return detected

# const
kDbPath = "../drosophilaDatabase_neop"

species_tree = "(TCAS,(((((BMOR,TRNI),(GMEL,HKAH)),((PPOL,PXUT),DPLE)),((CNAS,(AGAM,AALB)),(SLEB,(((DGRI,((DHYD,(DNAV,(DMOJ,DARI))),(DVIR,DNOV))),DBUS),(DWIL,(((DOBS,DGUA),((DPSE,DPER),DMIR)),((DBIP,DANA),((DSER,DKIK),((DELE,DRHO),((DTAK,(DSUZ,DBIA)),(DEUG,((DERE,DYAK),(DMEL,(DSIM,DSEC)))))))))))))),((MSAC,BTAB),(((BTER,AMEL),(OBIR,(MPHA,ACEP))),(CSOL,NVIT)))));"

# setup
handle = StringIO(species_tree)
tree = Phylo.read(handle, "newick")

c = get_db_cursor(kDbPath)

# temp
Phylo.draw_ascii(tree)

# db_temp = ["CG33687|NP_001027130.2|DMEL"]

# c.execute("SELECT id FROM processed_trees")
# db_genes = list(map(lambda gene: gene["id"], c.fetchall()))

# shork

db_genes = []

with open("./abs_res/Detection_failure_probabilities") as f:
  while True:
    l = f.readline().strip()

    if l.startswith("#"):
      continue

    if l is "":
      break

    sp = l.split("\t")

    if sp[0].find("DSUZ") == -1:
      continue
    
    db_genes.append(sp[0])

# a

# 1 gene can only be in 1 orthogroup (load ALL genes in file and compare len of list vs set)

print("total genes: {}".format(len(db_genes)))

curr = 0

tempa = []
tempb = []

with open("../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt") as file:
  while True:
    orthogroup = file.readline()

    if orthogroup == "":
      break

    genes = list(map(lambda gene: gene.strip(), orthogroup[(orthogroup.find(":") + 2):].split(" ")))
    species_uniq = list(set(map(lambda gene: gene[(gene.rfind("|") + 1):], genes)))

    for db_gene in db_genes:

      if db_gene not in genes:
        continue

      lowest_path = sys.maxsize

      # lazy and inefficient but it works so ye
      # well its not exactly that but there *is* a better way, its just not worth the effort
      for sp in species_uniq:

        db_gene_sp = db_gene[(db_gene.rfind("|") + 1):]

        mrca = tree.common_ancestor(sp, db_gene_sp)
        path = tree.get_path(mrca)

        lowest_path = len(path) if len(path) < lowest_path else lowest_path

      curr += 1

      if lowest_path == 0:
        continue

      # check if detected (absense)
      det_count = homolog_det_count(db_gene)

      if det_count == 0:
        print(db_gene)
        continue

      tempa.append(lowest_path)
      tempb.append(det_count)

      print("{}. {}: {} - {}".format(curr, db_gene, lowest_path, det_count))

from matplotlib import pyplot
from scipy.stats import pearsonr, spearmanr

pyplot.scatter(tempa, tempb)
pyplot.savefig("trg_abs.png")

pcor, _ = pearsonr(tempa, tempb)
scor, _ = spearmanr(tempa, tempb)

print("pearson: {}\nspearman: {}".format(pcor, scor))
