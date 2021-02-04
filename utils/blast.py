import sys
import gzip

from functools import lru_cache as cache

def process_blast_line(line, gz = False):
  line = line.decode("utf-8") if gz else line

  parts = line.strip().split("\t")

  return {
    "query": parts[0],
    "subject": parts[1],
    "identity": float(parts[2]),
    "alignment": int(parts[3]),
    "mismatches": int(parts[4]),
    "gaps": int(parts[5]),
    "qstart": int(parts[6]),
    "qend": int(parts[7]),
    "sstart": int(parts[8]),
    "send": int(parts[9]),
    "e": float(parts[10]),
    "score": float(parts[11]),
  }

@cache(maxsize = None)
def read_blast_file(blast_fp, gz = False):
  entries = []

  open_fn = gzip.open if gz else open

  with open_fn(blast_fp) as file:
    entries = list(map(lambda e: process_blast_line(e, gz), file.readlines()))

  return entries

def best_hits(blast_fp):
  blast = read_blast_file(blast_fp)

  # list of lowest e value hits for x query
  lowest_hits = []

  # (query, subject, eval)
  last_query = ""

  for entry in blast:
    query = entry["query"]
    subject = entry["subject"]
    e = entry["e"]

    # if we moved onto a new query gene, reset last hit
    if last_query != query:
      last_query = query
      lowest_hits.append((query, subject, e))
      continue

    if e < lowest_hits[-1][2]:
      lowest_hits[-1] = (query, subject, e)

  return lowest_hits

# given two blast files, returns a list of reciprocal best hits
def rec_best_hits(blast_fp1, blast_fp2):
  lowest_hits1 = best_hits(blast_fp1)
  lowest_hits2 = best_hits(blast_fp2)

  lowest_reciprocal = []

  for q1, s1, _ in lowest_hits1:
    for q2, s2, _ in lowest_hits2:
      if (q1, s1) == (s2, q2):
        lowest_reciprocal.append(q1)

  return lowest_reciprocal
