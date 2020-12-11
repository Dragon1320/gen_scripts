from os import listdir
from os.path import isfile, join

# const
blastFile = "../../blast"
orthogroupPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"
orthogroupSeqPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroup_Sequences"

species = [
  "AALB",
  "ACEP",
  "AGAM",
  "AMEL",
  "BMOR",
  "BTAB",
  "BTER",
  "CNAS",
  "CSOL",
  "DANA",
  "DARI",
  "DBIA",
  "DBIP",
  "DBUS",
  "DELE",
  "DERE",
  "DEUG",
  "DGRI",
  "DGUA",
  "DHYD",
  "DKIK",
  "DMEL",
  "DMIR",
  "DMOJ",
  "DNAV",
  "DNOV",
  "DOBS",
  "DPER",
  "DPLE",
  "DPSE",
  "DRHO",
  "DSER",
  "DSIM",
  "DSUZ",
  "DTAK",
  "DVIR",
  "DWIL",
  "DYAK",
  "GMEL",
  "HKAH",
  "MPHA",
  "MSAC",
  "NVIT",
  "OBIR",
  "PPOL",
  "PXUT",
  "SLEB",
  "TCAS",
  "TRNI",
  "DSEC",
]

# read a file as dict, given a delimiter
def read_dict(path, delim, include = lambda key, val: True):
  res = {}

  with open(path) as file:
    # may be helpful for larger files (think 1000s upon 1000s of lines)
    # rather than file.readlines()
    while True:
      line = file.readline().replace("\\", "/")
      
      if line == "":
        break

      if line.find(delim) == -1:
        continue

      delim_loc = line.find(delim)

      key = line[:delim_loc].strip()
      val = line[delim_loc + 1:].strip()

      if not include(key, val):
        continue

      res[key] = val

  return res

def find_orthologs(path, genes):
  res = { k: [] for k in genes }

  with open(path) as file:
    while True:
      line = file.readline().replace("\\", "/")

      if line == "":
        break

      orthogroup = line[(line.find(":") + 1):].strip().split(" ")

      for gene in genes:
        if gene in orthogroup:
          res[gene] = orthogroup

  return res

blast = read_dict(blastFile, "\t")

genes = list(set(blast.keys()))

orthogorups = find_orthologs(orthogroupPath, genes)

# all the DMEL genes that have orthologs in all other species 
ortho_genes = []

for gk in orthogorups:
  orthogroup_species_list = list(set(map(lambda gene: gene[(gene.rfind("|") + 1):], orthogorups[gk])))

  if len(species) != len(orthogroup_species_list):
    continue

  ortho_genes.append(gk)

busco = []
for k, v in blast.items():
  if k not in ortho_genes:
    continue

  bg = v[:(v.find("\t"))]

  busco.append(bg)

busco = list(set(busco))

print(len(genes), len(ortho_genes), len(busco))

# print(ortho_genes[0])

# genes = {}
# for bg in busco:
#   for k, v in blast.items():
#     bg = line[:(line.find("\t"))]

#     if genes[bg] is None:
#       genes[bg] = []

#     genes[]

# print(busco, len(busco))

# exit(-1)

# print(len(genes))









# print(orthologs[genes[0]])

# orthogroup_ids = []

# orthogroup_file = read_dict(orthogroupPath, ":")

# for k in orthogroup_file.keys():
#   for g in genes:
#     if g not in orthogroup_file[k]:
#       continue

#     orthogroup_ids.append(k)

# print(len(orthogroup_ids))

# print(orthogroup_ids)

# get all seq files
# orthogroupSeqFiles = [f for f in listdir(orthogroupSeqPath) if isfile(join(orthogroupSeqPath, f)) and f.split(".fa")[0] in orthogroup_ids]

# fasta_file_paths = []

# for f in listdir(orthogroupSeqPath):
#   if isfile(join(orthogroupSeqPath, f)) and f.split(".fa")[0] in orthogroup_ids:
#     fasta_file_paths.append(join(orthogroupSeqPath, f))
