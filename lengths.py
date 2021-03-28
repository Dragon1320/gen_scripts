from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

from utils.data import get_fastevol_singletons, get_singletons, get_duplicable
from utils.const import translation_fasta_dir, database_fp
from utils.fasta import read_fasta

# TODO: move to utils?
def round_e(num):
  return float("{:.6e}".format(num))

def get_gene_lengths(gene_ids, fasta_file):
  lengths = []

  for id in gene_ids:
    gene_seq = fasta_file.get(id)

    if gene_seq == None:
      # debug - print anything thats in the database but not in the fasta file
      print(id)
      continue

    lengths.append(len(gene_seq))

  return lengths

# dup vs sing, sing vs fsing, dup vs fsing
dsuz_fasta_fp = "{}/DSUZ.longest_only.faa".format(translation_fasta_dir)
out_plt_fp = "./out/plots/plt_lengths.png"

dsuz_fasta = read_fasta(dsuz_fasta_fp)

fsing = get_fastevol_singletons(database_fp)
sing = get_singletons(database_fp)
dup = get_duplicable(database_fp)

print(len(fsing), len(sing), len(dup))

fsing_len = get_gene_lengths(fsing, dsuz_fasta)
sing_len = get_gene_lengths(sing, dsuz_fasta)
dup_len = get_gene_lengths(dup, dsuz_fasta)

dup_v_sing = mannwhitneyu(dup_len, sing_len, alternative = "two-sided")
sing_v_fsing = mannwhitneyu(sing_len, fsing_len, alternative = "two-sided")
dup_v_fsing = mannwhitneyu(dup_len, fsing_len, alternative = "two-sided")

print("dup_v_sing: {}\n".format(dup_v_sing))
print("sing_v_fsing: {}\n".format(sing_v_fsing))
print("dup_v_fsing: {}\n".format(dup_v_fsing))

# plot
fig, ax = plt.subplots()

fig.set_size_inches(10, 10)

ax.set_ylim(top = 5000)
ax.set_xlim(left = 0.5, right = 4)
ax.set_ylabel("gene length", fontsize = 16, labelpad = 20)
ax.set_xlabel("duplicability", fontsize = 16, labelpad = 20)
ax.tick_params(axis = "both", which = "major", labelsize = 12)
# ax.set_xlabel("group", fontsize = 11, fontweight = "bold")
# ax.set_xticks(["dup", "sing", "fsing"])
ax.set_xticklabels(["duplicate", "singleton", "fast-evolving singleton"])

ax.hlines([4000, 4000, 4500], [1, 2.05, 1], [1.95, 3, 3], color = "black")
ax.text(1.5, 4100, "p={}".format(round_e(dup_v_sing[1])), horizontalalignment = "center", fontsize = 12)
ax.text(2.5, 4100, "p={}".format(round_e(sing_v_fsing[1])), horizontalalignment = "center", fontsize = 12)
ax.text(2, 4600, "p={}".format(round_e(dup_v_fsing[1])), horizontalalignment = "center", fontsize = 12)

bp = ax.boxplot([dup_len, sing_len, fsing_len], showmeans = True, meanline = True)

ax.legend([bp["medians"][0], bp["means"][0]], ["median", "mean"], fontsize = 12)
# fig.subplots_adjust(right = 1)

plt.savefig(out_plt_fp)
