protdistOutfileFp = "./outfile"
focalSpName = "DSUZ"
outfileFp = "./distances.txt"

sp = []
vals = []

with open(protdistOutfileFp) as f:
  h = int(f.readline().strip())
  r = False
  c = 0

  while True:
    l = f.readline()

    if len(l) == 0:
      break

    if l[0] != " ":
      sp.append(l.split(" ")[0])

    l = l.strip()

    if l.split(" ")[0] != focalSpName:
      if not r:
        continue
    else:
      r = True
      l = l[(len(focalSpName) + 1):]

    if r and c < h:
      cl = list(filter(lambda f: f != "", l.split(" ")))

      for v in cl:
        vals.append(v)
        c += 1

    if c >= h:
      r = False

# print(vals)

outfile = open(outfileFp, "w+")

for i in range(0, len(sp)):
  outfile.write("{}\t{}\n".format(sp[i], vals[i]))

outfile.close()
