from os import listdir
from os.path import isfile, join
from functools import reduce

from utils.db_init import get_db_cursor

# const
kDbPath = "../drosophilaDatabase_neop"
orthogroupPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"
speciesIdsPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/WorkingDirectory/SpeciesIDs.txt"
sequenceIdsPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/WorkingDirectory/SequenceIDs.txt"
blastFolderPath = "./blast_bit"

focalSpeciesName = "DSUZ.longest_only.faa"

outFile = "./out/bitscores.txt"

species = [
  "AALB.longest_only.faa",
  "ACEP.longest_only.faa",
  "AGAM.longest_only.faa",
  "AMEL.longest_only.faa",
  "BMOR.longest_only.faa",
  "BTAB.longest_only.faa",
  "BTER.longest_only.faa",
  "CNAS.longest_only.faa",
  "CSOL.longest_only.faa",
  "DANA.longest_only.faa",
  "DARI.longest_only.faa",
  "DBIA.longest_only.faa",
  "DBIP.longest_only.faa",
  "DBUS.longest_only.faa",
  "DELE.longest_only.faa",
  "DERE.longest_only.faa",
  "DEUG.longest_only.faa",
  "DGRI.longest_only.faa",
  "DGUA.longest_only.faa",
  "DHYD.longest_only.faa",
  "DKIK.longest_only.faa",
  "DMEL.longest_only.faa",
  "DMIR.longest_only.faa",
  "DMOJ.longest_only.faa",
  "DNAV.longest_only.faa",
  "DNOV.longest_only.faa",
  "DOBS.longest_only.faa",
  "DPER.longest_only.faa",
  "DPLE.longest_only.faa",
  "DPSE.longest_only.faa",
  "DRHO.longest_only.faa",
  "DSER.longest_only.faa",
  "DSIM.longest_only.faa",
  "DSUZ.longest_only.faa",
  "DTAK.longest_only.faa",
  "DVIR.longest_only.faa",
  "DWIL.longest_only.faa",
  "DYAK.longest_only.faa",
  "GMEL.longest_only.faa",
  "HKAH.longest_only.faa",
  "MPHA.longest_only.faa",
  "MSAC.longest_only.faa",
  "NVIT.longest_only.faa",
  "OBIR.longest_only.faa",
  "PPOL.longest_only.faa",
  "PXUT.longest_only.faa",
  "SLEB.longest_only.faa",
  "TCAS.longest_only.faa",
  "TRNI.longest_only.faa",
  "DSEC.longest_only.faa",
]

# setup
c = get_db_cursor(kDbPath)

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

# script

# load species name id mappings
# kv map, where k = id, v = name
speciesIdsToNames = read_dict(speciesIdsPath, ":")
# kv map, where k = name, v = id
speciesNamesToIds = { v: k for k, v in speciesIdsToNames.items() }
focalSpeciesId = speciesNamesToIds.get(focalSpeciesName)

# load gene name id mappings for genes of interest
# WHERE dup_status IS 'S'?
c.execute("SELECT id FROM processed_trees WHERE (dup_status IS 'S' OR dup_status IS 'D') AND excludedReason IS NULL")
dbGeneNames = list(map(lambda gene: gene["id"], c.fetchall()))

focalSpeciesOrthologs = find_orthologs(orthogroupPath, dbGeneNames)

# will work without this part (key.startswith(focalSpeciesId)), BUT will be incredibly slow
# kv map, where k = id, v = name
sequenceIdsToNames = read_dict(sequenceIdsPath, ":")
# ye so we cant unpack these like above with speciesIds cos python is too fragile to handle over 2k entries..
# kv map, where k = name, v = id
sequenceNamesToIds = dict(zip(sequenceIdsToNames.values(), sequenceIdsToNames.keys()))

# get all the blast files
blastFiles = [f for f in listdir(blastFolderPath) if isfile(join(blastFolderPath, f)) and f.startswith("Blast{}".format(focalSpeciesId))]

# ugh
results = dict(zip(dbGeneNames, list(map(lambda db_gene: { k: [] for k in species }, dbGeneNames))))

# debug
missingKeys = []

for blastFile in blastFiles:
  with open("{}/{}".format(blastFolderPath, blastFile)) as file:

    print("processing file: {}".format(blastFile))

    while True:
      line = file.readline().strip().split("\t")

      if line[0] == "":
        break

      #
      focalSpeciesEntry = line[0]
      focalSpeciesId = focalSpeciesEntry.split("_")[0]
      focalSpeciesGeneId = focalSpeciesEntry.split("_")[1]

      compSpeciesEntry = line[1]
      compSpeciesId = compSpeciesEntry.split("_")[0]
      compSpeciesGeneId = compSpeciesEntry.split("_")[1]

      bitscore = line[-1]

      # get focal species gene name for current line
      focalSpeciesGeneName = sequenceIdsToNames["{}_{}".format(focalSpeciesId, focalSpeciesGeneId)]

      # skip if gene name not in the db
      if focalSpeciesGeneName not in dbGeneNames:
        continue

      # get all ortholog ids
      focalSpeciesGeneOrthologs = focalSpeciesOrthologs[focalSpeciesGeneName]

      focalSpeciesGeneOrthologIds = []
      for geneName in focalSpeciesGeneOrthologs:
        if geneName in sequenceNamesToIds.keys():
          focalSpeciesGeneOrthologIds.append(sequenceNamesToIds[geneName])
        else:
          # print("missing key: {}".format(geneName))
          missingKeys.append(geneName)

      # only check line if one of orthologs
      if "{}_{}".format(compSpeciesId, compSpeciesGeneId) in focalSpeciesGeneOrthologIds:
        results[focalSpeciesGeneName][speciesIdsToNames[compSpeciesId]].append(bitscore)
        # print(speciesIdsToNames[compSpeciesId])
        # print("{} - {} - {}".format(focalSpeciesGeneName, speciesIdsToNames[compSpeciesId], bitscore))

for dbgene in results:
  for sp in results[dbgene]:
    if len(results[dbgene][sp]) == 0:
      results[dbgene][sp] = 0
    else:
      results[dbgene][sp] = str(float(reduce(lambda a, b: str(float(a) + float(b)), results[dbgene][sp], "0")) / len(results[dbgene][sp]))

outFile = open(outFile, "w+")
outFile.write("{}\t{}\n".format("None", "\t".join(list(map(lambda sp: sp[:(sp.find("."))], species)))))

for dbgene in results:
  outFile.write("{}\t{}\n".format(dbgene, "\t".join(results[dbgene].values())))

outFile.close()

missingKeys = list(set(missingKeys))
print(missingKeys, len(missingKeys))

# questions:
# missing gene ids
