import gzip
import math
import os

from threading import Thread, Lock

from utils import helpers, genes

# const
thread_count = int(os.environ.get("THREAD_COUNT", 1))

gene_ids_fp = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/WorkingDirectory/SequenceIDs.txt"
species_ids_fp = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/WorkingDirectory/SpeciesIDs.txt"
orthogroup_fp = "../Results_Jun25/WorkingDirectory/OrthoFinder/Results_Aug27/Orthogroups/Orthogroups.txt"
res_fp = "./out/para.txt"

fastevol = genes.get_fastevol()
single = genes.get_single()

# file
flock = Lock()
res_file = open(res_fp, "w+")

# cache
lock = Lock()
blast_cache = {}

cutoffs = [
  float("1.0e-10"),
  float("1.0e-9"),
  float("1.0e-8"),
  float("1.0e-7"),
  float("1.0e-6"),
]

def find_paralog_species(genes):
  for cutoff in cutoffs:

    for gene in genes:
      orthologs = helpers.find_gene_orthologs(orthogroup_fp, gene)

      hits = []

      for orth in orthologs:
        species = helpers.gene_to_sp(orth)

        gene_id = helpers.dict_file_vk(gene_ids_fp, orth)
        species_id = helpers.dict_file_vk(species_ids_fp, species)

        if gene_id == None or species_id == None:
          print("error finding gene/species if for: {}".format(orth))
          continue

        blast_fp = "./blast/Blast{}_{}.txt.gz".format(species_id, species_id)

        # cache
        blast_file = []

        with lock:
          blast_file = blast_cache.get(blast_fp)

          if blast_file == None:
            blast_file = helpers.read_blast(blast_fp)
            blast_cache[blast_fp] = blast_file

        det = False

        for g1, g2, e in blast_file:
          if gene_id != g1 and det == True:
            break

          if gene_id == g1:
            det = True
          else:
            continue

          if g1 == g2:
            continue

          if float(e) >= cutoff:
            continue

          hits.append(g1.split("_")[0])

      sp_hits = list(map(lambda e: helpers.dict_file_kv(species_ids_fp, e), set(hits)))

      with flock:
        res_file.write("{}\t{}\t{}\n".format(gene, cutoff, "\t".join(sp_hits)))

      print("{} ({}): {}".format(gene, cutoff, sp_hits))

if __name__ == "__main__":
  gene_list = fastevol[:10]

  print("total genes: {}".format(len(gene_list)))

  thread_pool = []
  genes_per_thread = math.ceil(len(gene_list) / thread_count)

  for i in range(0, thread_count):
    genes = gene_list[(i * genes_per_thread):((i + 1) * genes_per_thread)]

    t = Thread(target=find_paralog_species, args=[genes])

    t.start()
    thread_pool.append(t)

  for t in thread_pool:
    t.join()
