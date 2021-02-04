from functools import lru_cache as cache

@cache(maxsize = None)
def read_fasta(fasta_fp):
  entries = {}

  with open(fasta_fp) as file:
    while True:
      header = file.readline().strip()
      seq = file.readline().strip()

      if header == "" or header == None:
        break

      # validation
      if header[0] != ">":
        print("invalid fasta file: {}".format(fasta_fp))
        exit(1)

      entries[header[1:]] = seq

  return entries
