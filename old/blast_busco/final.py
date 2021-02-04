rbhFp = "./rbh.txt"
dmelFp = "./out.txt"
outFileFp = "./final.txt"

groups = {}

lastGene = ""

with open(rbhFp) as f:
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

dmelBusco = []

with open(dmelFp) as f:
  dmelBusco = list(map(lambda l: l.strip(), f.readlines()))

finalGroups = {}

for g in groups:
  if len(groups[g]) == 50 and g in dmelBusco:
    finalGroups[g] = groups[g]

outFile = open(outFileFp, "w+")

for g in finalGroups:
  outFile.write(">{}\n".format(g))

  for pair in finalGroups[g]:
    outFile.write("{}\n".format(pair))

outFile.close()
