#!/usr/bin/env bash
# Install EModelX + CryoAtom2 envs/weights, then queue FAST_FIRST GPU runs.
set -u
CONDA=/home/hemangini/miniconda3/etc/profile.d/conda.sh
ROOT=/home/hemangini/Documents/temp/mica_bench
LOG=$ROOT/logs/install_emodelx_cryoatom.log
source "$CONDA"
exec >>"$LOG" 2>&1
echo "START $(date)"

# --- EModelX env (modern torch; official 1.8/CUDA 11.1 will not run on this box) ---
if ! conda env list | awk '{print $1}' | grep -qx EModelX; then
  conda create -y -n EModelX python=3.10
fi
conda activate EModelX
python -m pip install -U pip gdown
python -m pip install torch --index-url https://download.pytorch.org/whl/cu128
python -m pip install biopython numpy pandas scipy mrcfile tqdm scikit-image networkx superpose3d
python -m pip install open3d atom3d || python -m pip install open3d || true

WEIGHTS=$ROOT/src/EModelX/models
mkdir -p "$WEIGHTS"
if ! ls "$WEIGHTS"/*.ckpt >/dev/null 2>&1; then
  cd "$WEIGHTS"
  python -m gdown --fuzzy 'https://drive.google.com/file/d/13BKzEBfL0uubYgcJJTGSZ-PkO9hXAB9X/view?usp=drive_link' -O emodelx_weights.zip || true
  if [[ -f emodelx_weights.zip ]]; then
    unzip -o emodelx_weights.zip
    rm -f emodelx_weights.zip
    # zip may nest files
    find . -name '*.ckpt' -exec mv -n {} "$WEIGHTS"/ \;
  fi
fi
ls -lh "$WEIGHTS"
conda deactivate

# --- CryoAtom2 ---
cd "$ROOT/src/CryoAtom"
if [[ ! -f CryoAtom2/checkpoint/CryoNet.pth ]]; then
  wget -c --no-check-certificate https://yanglab.qd.sdu.edu.cn/CryoAtom/download/checkpoints_v2.1.zip -O /tmp/checkpoints_v2.1.zip
  unzip -o /tmp/checkpoints_v2.1.zip -d "$ROOT/src/CryoAtom"
fi
if ! conda env list | awk '{print $1}' | grep -qx CryoAtom2; then
  conda env create -f linux.yml
fi
source "$CONDA"
conda activate CryoAtom2
cd "$ROOT/src/CryoAtom"
python setup.py install
conda deactivate

echo "INSTALL_OK $(date)"
# Queue runs: EModelX first (paper method), CryoAtom waits inside its script for GPU.
systemd-run --user --unit=mica-bench-emodelx.service --collect \
  --property=WorkingDirectory="$ROOT" --property=MemoryMax=48G --property=Restart=no \
  /bin/bash -c "source $CONDA; cd $ROOT && python3 9_run_emodelx.py >>$ROOT/logs/emodelx_detached.log 2>&1"
systemd-run --user --unit=mica-bench-cryoatom.service --collect \
  --property=WorkingDirectory="$ROOT" --property=MemoryMax=48G --property=Restart=no \
  /bin/bash -c "source $CONDA; cd $ROOT && python3 10_run_cryoatom.py >>$ROOT/logs/cryoatom_detached.log 2>&1"
echo "QUEUED $(date)"
