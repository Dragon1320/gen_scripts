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
# TODO: that plot thing zoe had in her presentation
fig, ax = plt.subplots()

ax.set_ylim(top = 5000)
ax.set_ylabel("gene length")
# ax.set_xlabel("group", fontsize = 11, fontweight = "bold")
# ax.set_xticks(["dup", "sing", "fsing"])
ax.set_xticklabels(["duplicable", "singleton", "fast-evolving singleton"])

ax.hlines([4000, 4000, 4500], [1, 2.05, 1], [1.95, 3, 3])
ax.text(1.5, 4100, "p={}".format(round_e(dup_v_sing[1])), horizontalalignment = "center", fontsize = 8)
ax.text(2.5, 4100, "p={}".format(round_e(sing_v_fsing[1])), horizontalalignment = "center", fontsize = 8)
ax.text(2, 4600, "p={}".format(round_e(dup_v_fsing[1])), horizontalalignment = "center", fontsize = 8)

ax.boxplot([dup_len, sing_len, fsing_len])

plt.savefig(out_plt_fp)
