#!/usr/bin/env bash
# Start benchmark jobs as user systemd units so Cursor/chat hangup does not kill them.
set -u
ROOT=/home/hemangini/Documents/temp/mica_bench
LOG=$ROOT/logs
CONDA=/home/hemangini/miniconda3/etc/profile.d/conda.sh
mkdir -p "$LOG" "$LOG/detached"

start_unit() {
  local name=$1
  local cmd=$2
  local unit="mica-bench-${name}.service"
  if systemctl --user is-active --quiet "$unit"; then
    echo "$name already active"
    return
  fi
  systemctl --user reset-failed "$unit" 2>/dev/null || true
  systemd-run --user \
    --unit="$unit" \
    --collect \
    --same-dir \
    --setenv=HOME="$HOME" \
    --property=WorkingDirectory="$ROOT" \
    --property=MemoryMax=180G \
    --property=Restart=no \
    /bin/bash -c "source $CONDA; $cmd >>$LOG/${name}_detached.log 2>&1"
  echo "started $name ($unit)"
}

start_unit cryo2struct "conda activate cryo2struct && cd $ROOT && python 4_run_cryo2struct.py"
start_unit mica "cd $ROOT && python3 6_run_mica.py"
start_unit map_to_model "cd $ROOT && python3 7_run_phenix_map_to_model.py"
echo DETACHED_OK
systemctl --user list-units 'mica-bench-*.service' --no-legend
