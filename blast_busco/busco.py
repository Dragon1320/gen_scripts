buscoDmelFp = "./blast_dmel_busco"
dmelBuscoFp = "./blast_busco_dmel"
outFileFp = "./out.txt"

buscoDmel = []
dmelBusco = []

with open(buscoDmelFp) as f:
  buscoDmel = list(map(lambda l: l.strip(), f.readlines()))

with open(dmelBuscoFp) as f:
  dmelBusco = list(map(lambda l: l.strip(), f.readlines()))

# { buscoGene: (dmelGene, eVal) }
buscoValues = {}
# { dmelGene: (buscoGene, eVal) }
dmelValues = {}

for line in buscoDmel:
  fields = line.split("\t")

  buscoGene = fields[0]
  dmelGene = fields[1]
  eVal = float(fields[-2])

  if buscoGene not in buscoValues or buscoValues[buscoGene][1] > eVal:
    buscoValues[buscoGene] = (dmelGene, eVal)

for line in dmelBusco:
  fields = line.split("\t")

  dmelGene = fields[0]
  buscoGene = fields[1]
  eVal = float(fields[-2])

  if dmelGene not in dmelValues or dmelValues[dmelGene][1] > eVal:
    dmelValues[dmelGene] = (buscoGene, eVal)

final = []
outFile = open(outFileFp, "w+")

for busco in buscoValues:
  dmel = buscoValues[busco][0]

  # print(dmelValues.keys())

  # print(busco, dmel)
  # exit(-1)

  if dmel not in dmelValues:
    print("*")
    continue

  if dmelValues[dmel][0] == busco:
    final.append(dmel)
    outFile.write("{}\n".format(dmel))

# for dmel in dmelValues:
#   busco = dmelValues[dmel][0]

#   # print(dmelValues.keys())

#   # print(busco, dmel)
#   # exit(-1)

#   if busco not in buscoValues:
#     print("*")
#     continue

#   if buscoValues[busco][0] == dmel:
#     final.append(busco)

# print(len(buscoValues), len(dmelValues), len(final))

# print(final)

outFile.close()

# print(len(final), len(list(set(final))))
