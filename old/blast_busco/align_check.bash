out_dir="./aligned"

for entry in "$out_dir"/*; do
  fgrep -o ">" "$out_dir/${entry##*/}" | wc -l
done
