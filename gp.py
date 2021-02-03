from gprofiler import GProfiler

#
import sqlite3

import numpy as np
from scipy.stats import mannwhitneyu

from utils.db_init import get_db_cursor
from utils.helpers import find_dmel_orthologs

# const
kDbPath = "../drosophilaDatabase_neop"
orthogroupPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"

#
c = get_db_cursor(kDbPath)

c.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")

rates = list(map(lambda a: a["proxy_rate"], c.fetchall()))

mean = np.mean(rates)
std = np.std(rates)

cutoff = mean + std * 2

c.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL AND proxy_rate > :cutoff", { "cutoff": cutoff })

fastevol_suz = list(map(lambda a: a["id"], c.fetchall()))

fastevol = list(map(lambda e: e[:(e.find("|"))], find_dmel_orthologs(orthogroupPath, fastevol_suz)))

# print("\n".join(fastevol))
# exit(1)

print(len(fastevol))

#
background = []
with open("../neopTranslations_forOrthofinder/DMEL.longest_only.faa") as f:
  while True:
    h = f.readline().strip()
    _ = f.readline()

    if h == "":
      break

    if not h.startswith(">"):
      print("incorrect header: {}".format(h))
      break

    g = h[1:]
    gid = g[:(g.find("|"))]

    background.append(gid)

#
single = []

with open("./blast_sin_01") as f:
  while True:
    l = f.readline().strip()

    if l == "":
      break

    single.append(l)

#
other = []

for bg in background:
  if bg not in single:
    other.append(bg)

#
gp = GProfiler(return_dataframe=True)

for (genes, fname) in [(single, "single"), (fastevol, "fastevol"), (other, "other")]:
  for i in range(0, 2):
    underr = False
    if i == 1:
      underr = True

    res = gp.profile(
      query = {
        fname: genes,
      },
      organism = "dmelanogaster",
      sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"],
      user_threshold = 0.05,
      all_results = False,
      ordered = False,
      no_evidences = True,
      combined = False,
      measure_underrepresentation = underr,
      no_iea = False,
      domain_scope = "annotated",
      numeric_namespace = "ENTREZGENE",
      significance_threshold_method = "fdr",
      # background = background,
    )

  res.to_csv("./out/gprofiler_{}_{}.csv".format(fname, "normal" if not underr else "underr"))
  print(res)

# res2 = gp.profile(
#   query = {
#     "fastevol": fastevol,
#   },
#   organism = "dmelanogaster",
#   sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"],
#   user_threshold = 0.05,
#   all_results = False,
#   ordered = False,
#   no_evidences = True,
#   combined = False,
#   measure_underrepresentation = False,
#   no_iea = False,
#   domain_scope = "annotated",
#   numeric_namespace = "ENTREZGENE",
#   significance_threshold_method = "fdr",
#   # background = background,
# )

# res2.to_csv("./out/gprofiler_fastevol.csv")
# print(res2)

# res3 = gp.profile(
#   query = {
#     "other": other,
#   },
#   organism = "dmelanogaster",
#   sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"],
#   user_threshold = 0.05,
#   all_results = False,
#   ordered = False,
#   no_evidences = True,
#   combined = False,
#   measure_underrepresentation = False,
#   no_iea = False,
#   domain_scope = "annotated",
#   numeric_namespace = "ENTREZGENE",
#   significance_threshold_method = "fdr",
#   # background = background,
# )

# res3.to_csv("./out/gprofiler_other.csv")
# print(res3)
