from gprofiler import GProfiler

#
background = []
with open("../neopTranslations_forOrthofinder/DMEL.longest_only.faa") as f:
  while True:
    h = f.readline().strip()
    _ = f.readline()

    if h == "":
      break

    if not h.startswith(">"):
      print("incorrect header: {}".format(h))
      break

    g = h[1:]
    gid = g[:(g.find("|"))]

    background.append(gid)

#
fastevol = []

with open("./blast_sin_0001") as f:
  while True:
    l = f.readline().strip()

    if l == "":
      break

    fastevol.append(l)

#
# other = []

# for bg in background:
#   if bg not in fastevol:
#     other.append(bg)

#
gp = GProfiler(return_dataframe=True)

res = gp.profile(
  query = {
    "fastevol": fastevol,
    # "other": other,
  },
  organism = "dmelanogaster",
  sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"],
  user_threshold = 0.05,
  all_results = False,
  ordered = False,
  no_evidences = True,
  combined = False,
  measure_underrepresentation = False,
  no_iea = False,
  domain_scope = "annotated",
  numeric_namespace = "ENTREZGENE",
  significance_threshold_method = "fdr",
  background = background,
)

res.to_csv("out.csv")

print(res)
