from matplotlib import pyplot

test_fp = "./para.txt"
sing_fp = "./para_sing.txt"
fastevo_fp = "./para_fastevo.txt"

para_num = []
cutoff = []

with open(test_fp) as f:
  while True:
    l = f.readline().strip()

    if l == "" or l == None:
      break

    content = l.split("\t")
    para = content[2:]

    para_num.append(len(para))
    cutoff.append(content[1])

pyplot.scatter(cutoff, para_num)
pyplot.savefig("out.png")
