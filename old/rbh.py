import sys
import gzip
from os import listdir
from os.path import isfile, join

# consts
blastDir = "./blast"
# DMEL
focalSpeciesId = 21
blastSpeciesCount = 50
# focalSpeciesGeneFile = "../neopTranslations_forOrthofinder/DMEL.longest_only.faa"
sequenceIdsPath = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/WorkingDirectory/SequenceIDs.txt"
blastFile = "../../blast"
outFile = "./out/rbh.txt"

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

geneIdToName = read_dict(sequenceIdsPath, ":")
geneNameToId = dict(zip(geneIdToName.values(), geneIdToName.keys()))

blast = read_dict(blastFile, "\t")
genes = list(set(blast.keys()))

# temp
# focalGenes = genes

# print(genes[0])

# focalGenes = []
# with open(focalSpeciesGeneFile) as f:
#   while True:
#     l = f.readline().strip()

#     if len(l) == 0:
#       break

#     if l[0] != ">":
#       continue

#     focalGenes.append(l[1:])

# blastFiles = []
# for f in listdir(blastDir):
#   d = join(blastDir, f)

#   if isfile(d) and str(focalSpeciesId) in f:
#     blastFiles.append(d)

i = 0

outFile = open(outFile, "w+")

for g in genes:
  i += 1

  jId = geneNameToId[g]

  spId = jId.split("_")[0]
  gId = jId.split("_")[1]

  print(">{} - {} / {}".format(g, i, len(genes)))
  outFile.write(">{}".format(g))

  for idx in range(0, blastSpeciesCount):
    first = sys.maxsize
    firstJId1 = ""
    firstJId2 = ""

    second = sys.maxsize
    secondJId1 = ""
    secondJId2 = ""

    other = ""

    pf = False

    fp = join(blastDir, "Blast{}_{}.txt.gz".format(str(focalSpeciesId), str(idx)))
    with gzip.open(fp) as f:
      while True:
        l = f.readline().decode("utf-8").strip()

        if len(l) == 0:
          break
        
        vals = l.split("\t")
        blastJId = vals[0]

        if jId != blastJId:
          if pf:
            break

          continue

        pf = True

        e = float(vals[-2])

        if first > e:
          otherBlastJId = vals[1]

          first = e
          firstJId1 = blastJId
          firstJId2 = otherBlastJId

    pf = False

    fp = join(blastDir, "Blast{}_{}.txt.gz".format(str(idx), str(focalSpeciesId)))
    with gzip.open(fp) as f:
      while True:
        l = f.readline().decode("utf-8").strip()

        if len(l) == 0:
          break
        
        vals = l.split("\t")
        blastJId = vals[0]

        if firstJId2 != blastJId:
          if pf:
            break

          continue

        pf = True
        
        e = float(vals[-2])

        if second > e:
          otherBlastJId = vals[1]

          second = e
          secondJId1 = blastJId
          secondJId2 = otherBlastJId

    if firstJId1 == secondJId2:
      # rbhIds.append(firstJId1, secondJId1)
      # temp
      # print(firstJId1, secondJId1)
      outFile.write("{} {}".format(firstJId1, secondJId1))

outFile.close()

# print(len(focalGenes))

# 21 vs x - find lowest
# x vs 21 - find lowest
# same?

# test = blastFiles[0]

# with gzip.open(test) as f:
#   l = f.readline().strip()

#   print(l)

# print(blastFiles)
