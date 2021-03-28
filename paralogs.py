import seaborn as sns
import pandas as pd
from numpy import arange, zeros
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

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

# # debug
# fsing_entries.append([{
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

# graphs
# data_fsing = []
# data_sing = []

x_size = 20
y_size = 10

df_fsing = pd.DataFrame(zeros((y_size, x_size)))
df_sing = pd.DataFrame(zeros((y_size, x_size)))

# data = [[0 for i in range(x_size)] for j in range(y_size)]

for i in reversed(range(0, x_size)):
  cutoff = 10 ** -i

  for entry in fsing_entries:
    paralogs = filter_paralogs(entry, cutoff)
    paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

    df_fsing[i][len(paralog_sp)] += 1

  for entry in sing_entries:
    paralogs = filter_paralogs(entry, cutoff)
    paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

    df_sing[i][len(paralog_sp)] += 1

vmin = min(df_fsing.values.min(), df_sing.values.min())
vmax = max(df_fsing.values.max(), df_sing.values.max())

# heatmap
plt.ylabel("num species")
plt.xlabel("cutoff")

# plt.ylim(bottom = -0.5, top = 9.5)
# plt.xlim(right = 19.5)

# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()

# plt.xscale("log")
# plt.scatter(data_x, data_y, color = colours)

# plt.imshow(data, cmap = "hot", interpolation = "nearest")

fig, axs = plt.subplots(ncols = 3, gridspec_kw = dict(width_ratios = [2, 2, 0.1]))

fig.set_size_inches(16, 7)

# axs[0].ylabel("non-self hits", fontsize = 16, labelpad = 20)
# axs[0].xlabel("e-value cutoff", fontsize = 16, labelpad = 20)

sns.heatmap(df_fsing, cmap = plt.cm.Blues, cbar = False, ax = axs[0], vmin = vmin)
sns.heatmap(df_sing, cmap = plt.cm.Blues, cbar = False, yticklabels = False, ax = axs[1], vmax = vmax)

fig.colorbar(axs[1].collections[0], cax = axs[2])

# axs[0].set_ylim(bottom = -0.5, top = 9.5)
# axs[0].set_xlim(right = 19.5)
axs[0].invert_xaxis()
axs[0].invert_yaxis()
axs[0].set_title("fast-evolving singletons", fontsize = 16, pad = 20)
axs[0].tick_params(axis = "both", which = "major", labelsize = 12)

ticks = []

for tick in axs[0].get_xticklabels():
  ticks.append("{}{}".format("" if tick.get_text() == "0" else "-", tick.get_text()))

axs[0].set_xticklabels(ticks)

axs[0].locator_params(axis = "x", nbins = 10)

# axs[1].set_ylim(bottom = -0.5, top = 9.5)
# axs[1].set_xlim(right = 19.5)
axs[1].invert_xaxis()
axs[1].invert_yaxis()
axs[1].set_title("other singletons", fontsize = 16, pad = 20)
axs[1].tick_params(axis = "both", which = "major", labelsize = 13)

ticks = []

for tick in axs[1].get_xticklabels():
  ticks.append("{}{}".format("" if tick.get_text() == "0" else "-", tick.get_text()))

axs[1].set_xticklabels(ticks)

axs[1].locator_params(axis = "x", nbins = 10)

axs[2].tick_params(axis = "both", which = "major", labelsize = 12)
axs[2].get_yaxis().set_major_locator(ticker.MaxNLocator(integer = True))

fig.add_subplot(111, frameon = False, label = "a")
plt.tick_params(labelcolor = "none", top = False, bottom = False, left = False, right = False)
plt.grid(False)
plt.xlabel("e-value cutoff (log scale)", fontsize = 16, labelpad = 20)
plt.ylabel("species count (non-self)", fontsize = 16, labelpad = 20)

fig.add_subplot(111, frameon = False, label = "bottom")
plt.tick_params(labelcolor = "none", top = False, bottom = False, left = False, right = False)
plt.grid(False)
plt.gca().yaxis.set_label_position("right")
plt.ylabel("hit count", fontsize = 16, labelpad = 50)

plt.savefig(out_plt_fp)

# reset plot
# plt.clf()

# old...

# f = open("./temp.txt", "w+")

# for entry in (fsing_entries + sing_entries):
#   paralogs = filter_paralogs(entry, 1e-10)

#   for paralog in paralogs:
#     # species, query, subject, og, e
#     f.write("{}\t{}\t{}\t{}\t{}\n".format(
#       kv_lookup(sp_ids_fp, paralog["query"].split("_")[0]),
#       kv_lookup(gene_ids_fp, paralog["query"]),
#       kv_lookup(gene_ids_fp, paralog["subject"]),
#       find_gene_orthogroup(database_fp, kv_lookup(gene_ids_fp, paralog["query"])),
#       paralog["e"],
#     ))

# f.close()

# exit(1)

# for i in [21, 16, 11, 6, 1]:
#   cutoff = 10 ** -i

#   for entry in fsing_entries:
#     print()

#     paralogs = filter_paralogs(entry, cutoff)
#     paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

#     data_fsing.append(len(paralog_sp))

#   for entry in sing_entries:
#     paralogs = filter_paralogs(entry, cutoff)
#     paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

#     data_sing.append(len(paralog_sp))

#   plt.ylim(top = 10)
#   plt.hist([data_fsing, data_sing], color = ["red", "blue"], bins = arange(6), edgecolor = "gray")

#   plt.savefig("./out/plots/plt_paralogs_{}.png".format(i))

# # old graphs
# data_x = []
# data_y = []
# colours = []

# for i in reversed(range(0, 20)):
#   cutoff = 10 ** -i

#   for entry in fsing_entries:
#     paralogs = filter_paralogs(entry, cutoff)
#     paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

#     data_x.append(cutoff)
#     data_y.append(len(paralog_sp))
#     colours.append("#ff000077")
#     # colours.append("#0000ff77")

  # for entry in sing_entries:
  #   paralogs = filter_paralogs(entry, cutoff)
  #   paralog_sp = list(set(map(lambda e: e["query"].split("_")[0], paralogs)))

  #   data_x.append(cutoff)
  #   data_y.append(len(paralog_sp))
  #   # colours.append("#ff000077")
  #   colours.append("#0000ff77")

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

# # plot
# plt.ylabel("num species")
# plt.xlabel("cutoff")

# plt.xscale("log")

# plt.scatter(data_x, data_y, color = colours)

# plt.savefig("{}_b.png".format(out_plt_fp))
