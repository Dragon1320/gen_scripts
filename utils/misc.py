from itertools import islice
from functools import lru_cache as cache

from utils.const import id_file_sep

def nth(iterable, n, default = None):
  return next(islice(iterable, n, None), default)

# ideally, would make specialised fns for this instead of caching the generic
# 'read_dict', but ive got enough ram and am not doing ridiculously complex calcs anyway
@cache(maxsize = None)
def read_dict(fp, sep):
  res = {}

  with open(fp) as f:
    while True:
      l = f.readline()
      
      if l == "" or l == None:
        break

      if l.find(sep) == -1:
        continue

      delim_loc = l.find(sep)

      k = l[:delim_loc].strip()
      v = l[delim_loc + 1:].strip()

      res[k] = v

  return res

# @cache(maxsize = None)
def dict_reverse(dict_to_rev):
  return dict(zip(dict_to_rev.values(), dict_to_rev.keys()))

def gene_to_sp(gene):
  return gene[(gene.rfind("|") + 1):]

def kv_lookup(table_fp, key):
  lookup_table = read_dict(table_fp, id_file_sep)

  return lookup_table.get(key)

def vk_lookup(table_fp, value):
  lookup_table = dict_reverse(read_dict(table_fp, id_file_sep))

  return lookup_table.get(value)
