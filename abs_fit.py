# NOTE: run the bash script with the same name
# this only exists cos ~~i need numpy~~ python is convenient

distance_file_fp = "./tmp/distances.txt"
normal_out_fp = "./tmp/out/normal"
includeonly_out_fp = "./tmp/out/includeonly"
bitscore_file = "Predicted_bitscores"
higherbound_bitscore_file = "Bitscore_99PI_higherbound_predictions"
lowerbound_bitscore_file = "Bitscore_99PI_lowerbound_predictions"

def parse_bitscore_entry(bitscore_entry):
  bitscore = None

  if "(Det)" in bitscore_entry:
    bitscore = bitscore_entry[:-5]
  elif "Det" in bitscore_entry:
    bitscore = bitscore_entry[(bitscore_entry.find(":") + 1):(bitscore_entry.find(")"))]
  else:
    bitscore = bitscore_entry

  return bitscore

def read_distance_file(fp):
  lines = []

  with open(fp) as file:
    lines = list(map(lambda e: e.strip().split("\t"), filter(lambda e: e != "", file.readlines())))

  lines.sort(key = lambda e: e[1])

  return lines

def read_bitscore_file(fp):
  species = None
  bitscores = None

  with open(fp) as file:
    while True:
      line = file.readline().strip()

      # eof
      if line == None or line == "":
        break

      # comment
      if line.startswith("#"):
        continue

      # header
      if species == None:
        species = line.split("\t")[1:]

        continue

      # distances (ignore everything after the first line)
      bitscores = list(map(parse_bitscore_entry, line.split("\t")[1:]))

      break

  merged = dict(zip(species, bitscores))

  return merged

# check if the predicted bitscore lies within the real higher/lower bounds
# all bitscores *should* be in the same order
bitscore_predicted= read_bitscore_file("{}/{}".format(includeonly_out_fp, bitscore_file))
bitscore_real_higher = read_bitscore_file("{}/{}".format(normal_out_fp, higherbound_bitscore_file))
bitscore_real_lower = read_bitscore_file("{}/{}".format(normal_out_fp, lowerbound_bitscore_file))

for sp, bitscore in bitscore_real.items():
  bitscore_predicted = float(bitscore)
  bitscore_real_higher = float(bitscore_real_higher[sp])
  bitscore_real_lower = float(bitscore_real_lower[sp])

  # print(bitscore, bitscore_higher, bitscore_lower)

  if (bitscore_real < bitscore_lower) or (bitscore_real > bitscore_higher):
    print("non-matching bitscore - bitscore: {}, higher: {}, lower: {}, species: {}".format(bitscore, bitscore_higher, bitscore_lower, sp))
