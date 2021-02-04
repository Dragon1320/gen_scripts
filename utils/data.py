from functools import lru_cache as cache

import numpy as np

from utils.db import connect

@cache(maxsize = None)
def get_all_genes(database_fp):
  cursor = connect(database_fp)

  # exclude all without data or with or with an excluded field
  cursor.execute("SELECT id FROM processed_trees WHERE excludedReason IS NULL AND proxy_rate IS NOT NULL")

  genes = list(map(lambda e: e["id"], cursor.fetchall()))

  return genes

@cache(maxsize = None)
def get_fastevol_singletons(database_fp):
  cursor = connect(database_fp)

  cursor.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")

  rates = list(map(lambda e: e["proxy_rate"], cursor.fetchall()))

  mean = np.mean(rates)
  std = np.std(rates)

  # hardcode for now, but its the most convenient, shouldnt need to change it either so its fiiine
  cutoff = mean + std * 2

  cursor.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL AND proxy_rate > :cutoff", { "cutoff": cutoff })

  fastevol = list(map(lambda e: e["id"], cursor.fetchall()))

  return fastevol

@cache(maxsize = None)
def get_singletons(database_fp):
  cursor = connect(database_fp)

  cursor.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")

  rates = list(map(lambda e: e["proxy_rate"], cursor.fetchall()))

  mean = np.mean(rates)
  std = np.std(rates)

  cutoff = mean + std * 2

  cursor.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL AND proxy_rate <= :cutoff", { "cutoff": cutoff })

  singletons = list(map(lambda e: e["id"], cursor.fetchall()))

  return singletons

@cache(maxsize = None)
def get_duplicable(database_fp):
  cursor = connect(database_fp)

  cursor.execute("SELECT id FROM processed_trees WHERE dup_status IS 'D' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")

  duplicable = list(map(lambda e: e["id"], cursor.fetchall()))

  return duplicable
