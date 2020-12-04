# const
blastFile = "../../blast"
orthogroupPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"

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

genes = list(read_dict(blastFile, "\t").keys());

orthogorups = find_orthologs(orthogroupPath, genes)

for gene in orthogorups:
  orthogroup = list(set(map(lambda gene: gene[(gene.rfind("|") + 1):], orthogorups[gene])))

  if len(species) != len(orthogroup):
    continue

  print(gene)

# print(orthologs[genes[0]])
