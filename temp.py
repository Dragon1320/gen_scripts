speices = {}

with open("./temp.txt") as f:
  for l in f.readlines():
    entries = l.strip().split("\t")

    sp = entries[0]

    if speices.get(sp[:4]) == None:
      speices[sp[:4]] = 0

    speices[sp[:4]] += 1

print(speices)
