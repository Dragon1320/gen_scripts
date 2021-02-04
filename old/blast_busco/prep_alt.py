from os.path import join

finalFp = "./final.txt"
geneIdFp = "./SequenceIDs.txt"
spIdFp = "./SpeciesIDs.txt"
seqsBaseFp = "./genes"
outBaseFp = "./groups"

# read a file as dict, given a delimiter
def read_dict(path, delim, include = lambda key, val: True):
  res = {}

  with open(path) as file:
    # may be helpful for larger files (think 1000s upon 1000s of lines)
    # rather than file.readlines()
    while True:
      line = file.readline()#.replace("\\", "/")
      
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

#
geneIdToName = read_dict(geneIdFp, ":")
geneNameToId = dict(zip(geneIdToName.values(), geneIdToName.keys()))

spIdToName = read_dict(spIdFp, ":")
spNameToId = dict(zip(spIdToName.values(), spIdToName.keys()))

#
groups = {}

lastGene = ""

with open(finalFp) as f:
  while True:
    l = f.readline()

    if len(l) == 0:
      break

    l = l.strip()

    if len(l) == 0:
      continue

    if l[0] == ">":
      lastGene = l[1:]
      groups[lastGene] = []

      continue

    groups[lastGene].append(l)

# ik this is probs gonna be slow af but its the easier way to do it
# if i was going for speed, id rewrite it in rust
for g in groups:
  print(g) # TODO: note progress

  for gg in groups[g]:
    target = gg.split(" ")[1]
    spName = spIdToName[target.split("_")[0]]
    gName = geneIdToName[target]
    spName = gName[(gName.rfind("|") + 1):]

    groupsFile = open(join(outBaseFp, "{}.faa".format(spName)), "a+")

    with open(join(seqsBaseFp, "{}.longest_only.faa".format(spName))) as f:
      while True:
        l = f.readline().strip()

        if len(l) == 0:
          break

        if l[0] == ">" and l[1:] == gName:
          nl = f.readline().strip()

          groupsFile.write("{}\n".format(l))
          groupsFile.write("{}\n".format(nl))

          break

  groupsFile.close()
