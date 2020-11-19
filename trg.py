from io import StringIO
from functools import reduce

from Bio import Phylo

from utils.db_init import get_db_cursor

# const
kDbPath = "../drosophilaDatabase_neop"

# setup
treedata = "((Dmel, (Dsim, Dsec)), (Dere, Dyak))"
handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")

c = get_db_cursor(kDbPath)

# temp
Phylo.draw_ascii(tree)

# calculate how taxonomically restricted duplicable genes are
# 0 - mrca is root
# +1 for every level down
# species in brackets '()', only 1 means terminal node
c.execute("SELECT id, dupInSp FROM processed_trees WHERE dup_status IS 'D' AND excludedReason IS NULL")

a = []

for row in c.fetchall():
  sp = row["dupInSp"].split(",")

  mrca = tree.common_ancestor(sp)
  path = tree.get_path(mrca)

  print("{} - {} - ({})".format(row["id"], len(path), ",".join(sp)))
