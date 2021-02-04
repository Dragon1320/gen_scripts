import gzip

# utility
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

def read_blast(fp):
  res = []

  with gzip.open(fp) as f:
    while True:
      l = f.readline().decode("utf-8").strip()

      if l == "" or l == None:
        break
      
      fields = l.split("\t")

      res.append((fields[0], fields[1], fields[-2]))

  return res

def gene_to_sp(gene):
  return "{}.longest_only.faa".format(gene[(gene.rfind("|") + 1):])

# species/gene ids
dict_file_sep = ":"

dict_cache = {}

def dict_reverse(dict_to_rev):
  return dict(zip(dict_to_rev.values(), dict_to_rev.keys()))

def dict_file_kv(dict_file_fp, key):
  global dict_file_sep
  global dict_cache

  cached_dict_file = dict_cache.get(dict_file_fp)

  if cached_dict_file != None:
    return cached_dict_file.get(key)

  dict_file_new = read_dict(dict_file_fp, dict_file_sep)
  dict_cache[dict_file_fp] = dict_file_new

  return dict_file_new.get(key)

def dict_file_vk(dict_file_fp, value):
  global dict_file_sep
  global dict_cache

  cached_dict_file = dict_cache.get(dict_file_fp)

  if cached_dict_file != None:
    return dict_reverse(cached_dict_file).get(value)

  dict_file_new = read_dict(dict_file_fp, dict_file_sep)
  dict_cache[dict_file_fp] = dict_file_new

  return dict_reverse(dict_file_new).get(value)

# orthologs
orth_file_sep = ":"

orth_cache = {}

def find_gene_orthologs(ortholog_file_fp, gene):
  global orth_file_sep
  global orth_cache

  ortholog_file = orth_cache.get(ortholog_file_fp)

  if ortholog_file == None:
    ortholog_file = read_dict(ortholog_file_fp, orth_file_sep)
    orth_cache[ortholog_file_fp] = ortholog_file

  for og_genes in ortholog_file.values():
    og_gene_list = og_genes.split(" ")

    if gene in og_gene_list:
      return og_gene_list

  return []

# old

def find_dmel_orthologs(path, genes):
  res = []

  with open(path) as file:
    while True:
      line = file.readline()

      if line == "":
        break

      orthogroup = line[(line.find(":") + 1):].strip().split(" ")

      for gene in genes:
        if gene in orthogroup:
          for orth_gene in orthogroup:
            if orth_gene[-4:] == "DMEL":
              res.append(orth_gene)

  return res
