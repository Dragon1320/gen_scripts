from gprofiler import GProfiler

from utils.data import get_fastevol_singletons, get_singletons, get_duplicable
from utils.const import database_fp, orthogroups_fp
from utils.orth import find_species_orthologs

def strip_gene_postfix(gene):
  stripped = gene[:(gene.find("|"))]

  return stripped

gp_species = "DMEL"
gp_out_dir = "./out/gprofiler"

fsing = get_fastevol_singletons(database_fp)
sing = get_singletons(database_fp)
dup = get_duplicable(database_fp)

fsing_dmel = find_species_orthologs(orthogroups_fp, fsing, gp_species)
sing_dmel = find_species_orthologs(orthogroups_fp, sing, gp_species)
dup_dmel = find_species_orthologs(orthogroups_fp, dup, gp_species)

fsing_dmel_stripped = list(map(lambda e: strip_gene_postfix(e), fsing_dmel))
sing_dmel_stripped = list(map(lambda e: strip_gene_postfix(e), sing_dmel))
dup_dmel_stripped = list(map(lambda e: strip_gene_postfix(e), dup_dmel))

gp = GProfiler(return_dataframe=True)

datasets = [
  ("fsing", fsing_dmel_stripped, False),
  ("fsing_ns", fsing_dmel_stripped, True),
  ("sing", sing_dmel_stripped, False),
  ("dup", dup_dmel_stripped, False),
]

for name, data, ns in datasets:
  for i in range(0, 2):
    underr = False

    if i == 1:
      underr = True
    
    res = gp.profile(
      query = {
        name: data,
      },
      organism = "dmelanogaster",
      sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"],
      user_threshold = 0.05,
      all_results = ns,
      ordered = False,
      no_evidences = True,
      combined = False,
      measure_underrepresentation = underr,
      no_iea = False,
      domain_scope = "annotated",
      numeric_namespace = "ENTREZGENE",
      significance_threshold_method = "fdr",
      # background = background,
    )

    res.to_csv("{}/gp_{}_{}.csv".format(gp_out_dir, name, "normal" if not underr else "underrepresented"))
    print(res)
