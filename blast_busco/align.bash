search_dir="./groups"
out_dir="./aligned"

for entry in "$search_dir"/*; do
  echo "${entry##*/}"
  muscle -in "$search_dir/${entry##*/}" -out "$out_dir/${entry##*/}"
done
