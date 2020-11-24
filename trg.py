import sys
from io import StringIO
from functools import reduce
from enum import Enum

from Bio import Phylo

from utils.db_init import get_db_cursor

# const
kDbPath = "../drosophilaDatabase_neop"

species_tree = "(TCAS,(((((BMOR,TRNI),(GMEL,HKAH)),((PPOL,PXUT),DPLE)),((CNAS,(AGAM,AALB)),(SLEB,(((DGRI,((DHYD,(DNAV,(DMOJ,DARI))),(DVIR,DNOV))),DBUS),(DWIL,(((DOBS,DGUA),((DPSE,DPER),DMIR)),((DBIP,DANA),((DSER,DKIK),((DELE,DRHO),((DTAK,(DSUZ,DBIA)),(DEUG,((DERE,DYAK),(DMEL,(DSIM,DSEC)))))))))))))),((MSAC,BTAB),(((BTER,AMEL),(OBIR,(MPHA,ACEP))),(CSOL,NVIT)))));"

# setup
handle = StringIO(species_tree)
tree = Phylo.read(handle, "newick")

c = get_db_cursor(kDbPath)

# temp
Phylo.draw_ascii(tree)

db_temp = ["CG33687|NP_001027130.2|DMEL"]

c.execute("SELECT id FROM processed_trees")
db_genes = list(map(lambda gene: gene["id"], c.fetchall()))

# 1 gene can only be in 1 orthogroup (load ALL genes in file and compare len of list vs set)

print("total genes: {}".format(len(db_genes)))

curr = 0

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

      print("{}. {}: {}".format(curr, db_gene, lowest_path))
