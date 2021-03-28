from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

from utils.db import connect
from utils.const import database_fp

out_plt_fp = "./out/plots/plt_rates.png"

# TODO: move to utils?
def round_e(num):
  return float("{:.6e}".format(num))

cursor = connect(database_fp)

cursor.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'S' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")
sing_rates = list(map(lambda e: e["proxy_rate"], cursor.fetchall()))

cursor.execute("SELECT proxy_rate FROM processed_trees WHERE dup_status IS 'D' AND excludedReason IS NULL AND proxy_rate IS NOT NULL")
dup_rates = list(map(lambda e: e["proxy_rate"], cursor.fetchall()))

sing_vs_dup = mannwhitneyu(sing_rates, dup_rates, alternative = "two-sided")

print("sing_vs_dup: {}\n".format(sing_vs_dup))

fig, ax = plt.subplots()

fig.set_size_inches(10, 10)

ax.set_ylim(top = 1.8)
ax.set_xlim(left = 0.5, right = 3)
ax.set_ylabel("dN/dS", fontsize = 16, labelpad = 20)
ax.set_xlabel("duplicability", fontsize = 16, labelpad = 20)
ax.tick_params(axis = "both", which = "major", labelsize = 12)
# ax.set_xlabel("group", fontsize = 11, fontweight = "bold")
# ax.set_xticks(["dup", "sing", "fsing"])
ax.set_xticklabels(["singleton", "duplicate"])

ax.hlines([1.6], [1], [2], color = "black")
ax.text(1.5, 1.63, "p={}".format(round_e(sing_vs_dup[1])), horizontalalignment = "center", fontsize = 12)

bp = ax.boxplot([sing_rates, dup_rates], showmeans = True, meanline = True)

ax.legend([bp["medians"][0], bp["means"][0]], ["median", "mean"], fontsize = 12)
# fig.subplots_adjust(right = 1)

plt.savefig(out_plt_fp)

