from os import listdir
from os.path import isfile, join

from utils.db_init import get_db_cursor

# const
kDbPath = "../drosophilaDatabase_neop"
outFile = "./out/lengths.txt"
focalSpeciesFastaFile = "../neopTranslations_forOrthofinder/DSUZ.longest_only.faa"

# setup
c = get_db_cursor(kDbPath)

# script
c.execute("SELECT id FROM processed_trees WHERE dup_status IS 'S'")

focalSpeciesGenes = list(map(lambda entry: entry["id"], c.fetchall()))

outFile = open(outFile, "w+")
outFile.write("{}\t{}\n".format("GeneID", "GeneLength"))

with open(focalSpeciesFastaFile) as file:
  while True:
    (header, gene) = (file.readline().strip(), file.readline().strip())

    if header == "":
      break

    # header validation
    if header[:1] != ">":
      raise Exception("expected fasta header '>', instead found {}".format(header[:1]))

    for focalSpeciesGene in focalSpeciesGenes:
      if header[1:] != focalSpeciesGene:
        continue

      outFile.write("{}\t{}\n".format(focalSpeciesGene, len(gene)))


outFile.close()
