from os import listdir
from os.path import isfile, join

alignedFp = "./aligned"
outfileFp = "../../out/infile"

alignedFiles = []

for f in listdir(alignedFp):
  if isfile(join(alignedFp, f)):
    alignedFiles.append(f)

seqs = {}

for fn in alignedFiles:
  fp = join(alignedFp, fn)

  with open(fp) as f:
    currSp = ""

    while True:
      l = f.readline().strip()

      if len(l) == 0:
        break

      if l[0] == ">":
        g = l[1:]
        sp = g[(g.rfind("|") + 1):]

        currSp = sp

        if sp not in seqs:
          seqs[sp] = []

        continue

      # sanity (something im lacking rn) check
      if l[0] != ">" and currSp == "":
        print("lines out of order")
        exit(-1)

      if currSp not in seqs:
        seqs[currSp] = []

      seqs[currSp].append(l)

seqs = dict(zip(seqs.keys(), list(map(lambda s: "".join(s), seqs.values()))))

count = len(seqs.keys())
lengths = len(list(seqs.values())[0])

outfile = open(outfileFp, "w+")

outfile.write("{} {}\n".format(count, lengths))

for s in seqs:
  outfile.write("{} {}\n".format(s.ljust(10), seqs[s]))

outfile.close()

print(count, lengths)
