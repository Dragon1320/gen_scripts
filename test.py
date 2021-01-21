hits = {}

with open("../dmel_self") as f:
  # pos = None

  while True:
    l = f.readline().strip()

    if l is "":
      break

    sp = l.split("\t")

    if hits.get(sp[0]) == None:
      hits[sp[0]] = []

    hits[sp[0]].append(sp[0])

single = list(filter(lambda e: len(e) == 1, list(hits.values())))
flat = [item for sublist in single for item in sublist]

print(len(flat))

with open("./blast_sin_0001", "w+") as f:
  for g in flat:
    gid = g[:(g.find("|"))]

    f.write("{}\n".format(gid))

# 100% of homolog detection in all species (true)
# def homolog_det(gene):
#   with open("./abs_res/Detection_failure_probabilities") as f:
#     while True:
#       l = f.readline().strip()

#       if l.startswith("#"):
#         continue

#       if l is "":
#         break

#       sp = l.split("\t")

#       if sp[0] != gene:
#         continue

#       for ip in range(1, len(sp)):
#         if sp[ip] == "Detected" or sp[ip] == "0.0":
#           continue

#         return False

#       return True

#   return False

# res = homolog_det("LOC108005253|XP_016923935.1|DSUZ")
# print(res)

# def homolog_det_count(gene):
#   with open("./abs_res/Detection_failure_probabilities") as f:
#     detected = 0

#     while True:
#       l = f.readline().strip()

#       if l.startswith("#"):
#         continue

#       if l is "":
#         break

#       sp = l.split("\t")

#       if sp[0] != gene:
#         continue

#       for ip in range(1, len(sp)):
#         if sp[ip] == "Detected" or sp[ip] == "0.0":
#           detected += 1

#       return detected

#   return detected

# res = homolog_det_count("LOC108013288|XP_016934548.1|DSUZ")
# print(res)
