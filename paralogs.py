from matplotlib import pyplot as plt

from utils.misc import kv_lookup, vk_lookup, gene_to_sp
from utils.const import blast_file_dir, ingroup_sp, orthogroups_fp, gene_ids_fp, sp_ids_fp, database_fp
from utils.orth import find_gene_orthologs, find_gene_orthogroup
from utils.blast import read_blast_file
from utils.data import get_fastevol_singletons, get_singletons

out_plt_fp = "./out/plots/plt_paralogs.png"

# paralong if: has non self hits at some e value cutoff
def find_paralogs(gene):
  orthologs = find_gene_orthologs(database_fp, gene)
  ingroup_orthologs = list(filter(lambda e: gene_to_sp(e) in ingroup_sp, orthologs))

  paralogs = []
  # TODO: this value is incorrect? i think?
  paralog_sp = 0

  for ortholog in ingroup_orthologs:
    species = gene_to_sp(ortholog)
    sp_id = vk_lookup(sp_ids_fp, "{}.longest_only.faa".format(species))
    gene_id = vk_lookup(gene_ids_fp, ortholog)

    blast_file_fp = "{}/Blast{}_{}.txt.gz".format(blast_file_dir, sp_id, sp_id)

    blast_file = read_blast_file(blast_file_fp, gz = True)

    non_self_hits = list(filter(lambda e: e["query"] == gene_id and e["query"] != e["subject"], blast_file))

    paralogs += non_self_hits

    if len(non_self_hits) > 0:
      paralog_sp += 1

  return (paralogs, paralog_sp)

def filter_paralogs(blast_entries, cutoff):
  threshold_hits = list(filter(lambda e: e["e"] < cutoff, blast_entries))

  return threshold_hits

# ik i can make a utility fn for this, i will at some point (probably)
fsing_entries = []
sing_entries = []

# debug
# entries.append([{
#   "query": "0_1",
#   "e": 1e-40,
# },{
#   "query": "0_2",
#   "e": 1e-10,
# },{
#   "query": "1_1",
#   "e": 1e-30,
# }])

# fastevo
fsing = get_fastevol_singletons(database_fp)
print(len(fsing))

for gene in fsing:
  (paralogs, paralog_sp) = find_paralogs(gene)

  if paralog_sp > 0:
    fsing_entries.append(paralogs)

  print(paralogs, paralog_sp)

# exit(1)

# test with 100 first (before commiting probs an hour to run this)
sing = get_singletons(database_fp)[:100]
print(len(sing))

for gene in sing:
  (paralogs, paralog_sp) = find_paralogs(gene)

  if paralog_sp > 0:
    sing_entries.append(paralogs)

  print(paralogs, paralog_sp)

data_x = []
data_y = []
colours = []

for i in reversed(range(0, 20)):
  cutoff = 10 ** -i

  for entry in fsing_entries:
    paralogs = filter_paralogs(entry, cutoff)
    paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

    data_x.append(cutoff)
    data_y.append(len(paralog_sp))
    colours.append("#ff000077")
    # colours.append("#0000ff77")

  for entry in sing_entries:
    paralogs = filter_paralogs(entry, cutoff)
    paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

    data_x.append(cutoff)
    data_y.append(len(paralog_sp))
    # colours.append("#ff000077")
    colours.append("#0000ff77")

# weird_ones = []

# cutoff = 1e-10

# for entry in fsing_entries:
#   for hit in entry:
#     if hit["e"] < cutoff:
#       weird_ones.append(hit)

# f = open("./out/paralogs.txt", "w+")

# species, query, subject, orthogroup (query), orthogroup (subject), e-val
# for w in weird_ones:
#   species = kv_lookup(sp_ids_fp, w["query"].split("_")[0])
#   query = kv_lookup(gene_ids_fp, w["query"])
#   subject = kv_lookup(gene_ids_fp, w["subject"])
#   orth_q = find_gene_orthogroup(orthogroups_fp, query)
#   orth_s = find_gene_orthogroup(orthogroups_fp, subject)
#   e = w["e"]

#   f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(species, query, subject, orth_q, orth_s, e))

# f.close()

# plot
plt.ylabel("num species")
plt.xlabel("cutoff")

plt.xscale("log")

plt.scatter(data_x, data_y, color = colours)

plt.savefig(out_plt_fp)

print(data_x, data_y)
