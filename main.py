import sqlite3

import numpy as np
from scipy.stats import mannwhitneyu

# utils and setup
def dict_factory(cursor, row):
  d = {}
  for idx, col in enumerate(cursor.description):
    d[col[0]] = row[idx]
  return d

conn = sqlite3.connect("../drosophilaDatabase_neop")
conn.row_factory = dict_factory
c = conn.cursor()

# calculate cutoff
c.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")

rates = list(map(lambda a: a["proxy_rate"], c.fetchall()))

mean = np.mean(rates)
std = np.std(rates)

cutoff = mean + std * 2

# get all ids over the cutoff
c.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL AND proxy_rate > :cutoff", { "cutoff": cutoff })

ids = list(map(lambda a: a["id"], c.fetchall()))

# get all other ids (singletons below cutoff and duplicables)
ids_other = []

c.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL AND proxy_rate <= :cutoff", { "cutoff": cutoff })

ids_other += list(map(lambda a: a["id"], c.fetchall()))

c.execute("SELECT id FROM processed_trees WHERE dup_status IS 'D' AND excludedReason IS NULL AND proxy_rate IS NOT NULL", { "cutoff": cutoff })

ids_other += list(map(lambda a: a["id"], c.fetchall()))

# look for ids and calc lengths
lengths = []
lengths_other = []

with open("../neopTranslations_forOrthofinder/DSUZ.longest_only.faa") as file:
  while True:
    (header, gene) = (file.readline().strip(), file.readline().strip())

    if header == None or header == "":
      break

    # header validation
    if header[:1] != ">":
      raise Exception("expected fasta header '>', instead found {}".format(header[:1]))

    header = header[1:]

    # check if file entry matches db id
    if header in ids:
      lengths.append(len(gene))

    if header in ids_other:
      lengths_other.append(len(gene))

# validation - check if all genes were found
if len(ids) != len(lengths):
  raise Exception("expected {} db genes, instead found {}".format(len(ids), len(lengths)))

if len(ids_other) != len(lengths_other):
  raise Exception("expected {} db genes (other), instead found {}".format(len(ids_other), len(lengths_other)))

res = mannwhitneyu(lengths, lengths_other, alternative="two-sided")

print(res)
