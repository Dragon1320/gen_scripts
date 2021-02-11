from functools import lru_cache as cache, reduce

from utils.misc import read_dict, gene_to_sp
from utils.const import orth_file_sep
from utils.db import connect

# @cache(maxsize = None)
# def read_orthogroup_file(orth_fp):
#   orth_file = read_dict(orth_fp, orth_file_sep)

#   # the read_dict fn reads the genes list as a long string,
#   # turn it into a list for convenience
#   values = []
#   values_str = orth_file.values()
#   for v in values_str:
#     values.append(v.split(" "))

#   return dict(zip(orth_file.keys(), values))

@cache(maxsize = None)
def read_orthogroup_file(database_fp):
  cursor = connect(database_fp)

  db_orth = list(cursor.execute("SELECT * FROM finalOrthogroups"))

  keys = list(map(lambda e: e["group_id"], db_orth))
  values = list(map(lambda e: e["assignedGenes"].split(","), db_orth))

  return dict(zip(keys, values))

def find_gene_orthologs(database_fp, gene):
  orth_file = read_orthogroup_file(database_fp)

  for og_gene in orth_file.values():
    if gene in og_gene:
      return og_gene

  return []

def find_species_orthologs(database_fp, genes, species):
  sp_orth = []

  for gene in genes:
    orth = find_gene_orthologs(database_fp, gene)

    sp_orth += list(filter(lambda e: gene_to_sp(e) == species, orth))

  return sp_orth


def find_gene_orthogroup(database_fp, gene):
  orth_file = read_orthogroup_file(database_fp)

  for og in orth_file.keys():
    if gene in orth_file[og]:
      return og
