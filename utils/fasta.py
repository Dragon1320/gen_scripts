from functools import lru_cache as cache

@cache(maxsize = None)
def read_fasta(fasta_fp):
  entries = {}

  header = ""
  seq = ""

  with open(fasta_fp) as file:
    while True:
      line = file.readline().strip()

      if line == "" or line == None:
        if len(header) > 0 and len(seq) > 0:
          entries[header] = seq

          header = ""
          seq = ""

        break

      if line[0] == ">":
        if len(header) > 0 and len(seq) > 0:
          entries[header] = seq

          header = ""
          seq = ""

        header = line[1:]
        continue

      seq += line

      # header = file.readline().strip()
      # seq = file.readline().strip()

      # if header == "" or header == None:
      #   break

      # # validation
      # if header[0] != ">":
      #   print("invalid fasta file: {}".format(fasta_fp))
      #   exit(1)

      # entries[header[1:]] = seq

  return entries
