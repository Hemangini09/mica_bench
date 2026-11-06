#!/usr/bin/env bash
# Official MICA path from the README (MICA_pipeline.sh): process AF3, dock, then
# run.py with the default AF3_struct protocol. Starts only after EModelX exits.
# 9TBI is not included.
set -u
ROOT=/home/hemangini/Documents/temp/mica_bench
MICA=$ROOT/src/MICA
CONDA=/home/hemangini/miniconda3/etc/profile.d/conda.sh
PHENIX=$ROOT/phenix/phenix-2.2.1-6174/phenix_env.sh
PULCHRA=$MICA/modules/pulchra304/src/pulchra
LOG=$ROOT/logs/mica_official.log

exec >>"$LOG" 2>&1
echo "==== official MICA waiter $(date) ===="

while systemctl --user is-active --quiet mica-bench-emodelx.service; do
  echo "waiting for EModelX $(date)"
  sleep 30
done
if pgrep -f 'run.py --protocol=temp_flex' >/dev/null; then
  echo "EModelX process still alive, waiting"
  while pgrep -f 'run.py --protocol=temp_flex' >/dev/null; do
    sleep 30
  done
fi

echo "EModelX done $(date). Starting official MICA pipeline for 9XON 9W22 9WSD"
source "$CONDA"
conda activate MICA
cd "$MICA"

# pdb:contour:resolution:emdb  (contour from EMDB author recommendation)
maps=(
  "9XON:0.0409:2.5:67084"
  "9W22:7.0:3.12:65563"
  "9WSD:0.15:4.0:66197"
)

failed=()
for spec in "${maps[@]}"; do
  IFS=: read -r pdb contour res emdb <<<"$spec"
  work=$ROOT/work/$pdb
  out=$ROOT/models/mica/$pdb
  mkdir -p "$out"
  echo "MICA official $pdb contour=$contour res=$res"
  before=$(date +%s)
  if ./MICA_pipeline.sh \
      -f "$work/$pdb.fasta" \
      -a "$work/AF3_results" \
      -m "$work/emd_${emdb}.mrc" \
      -c "$contour" \
      -r "$res" \
      -p "$PULCHRA" \
      -x "$PHENIX" \
      -d cuda; then
    echo "pipeline ok $pdb"
  else
    echo "pipeline FAIL $pdb"
    failed+=("$pdb")
  fi
  find "$MICA/output" "$work" -name '*all_atom_model.pdb' -newermt "@$before" -print -exec cp -a {} "$out/" \;
done
echo "mica official failed: ${failed[*]:-none}"
