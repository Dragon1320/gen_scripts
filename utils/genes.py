import numpy as np

from utils.helpers import find_dmel_orthologs
from utils.db_init import get_db_cursor

def get_fastevol():
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

  fastevol = list(map(lambda a: a["id"], c.fetchall()))

  # fastevol = list(map(lambda e: e[:(e.find("|"))], find_dmel_orthologs(orthogroupPath, fastevol_suz)))

  return fastevol

def get_single():
  single = []

  with open("./blast_sin_dsuz_01") as f:
    while True:
      l = f.readline().strip()

      if l == "":
        break

      single.append(l)

  return single

def get_background():
  background = []

  with open("../../neopTranslations_forOrthofinder/DMEL.longest_only.faa") as f:
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

  return background
